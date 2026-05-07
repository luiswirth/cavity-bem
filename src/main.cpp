#include <Bembel/AnsatzSpace>
#include <Bembel/Geometry>
#include <Bembel/H2Matrix>
#include <Bembel/IO>
#include <Bembel/LinearForm>
#include <Bembel/Maxwell>
#include <Eigen/Dense>
#include <unsupported/Eigen/IterativeSolvers>

#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

static constexpr double PI = 3.14159265358979323846;

using namespace std::complex_literals;
using complex = std::complex<double>;

complex green_scalar(double r, double k) {
  return std::exp(1i * k * r) / (4 * PI * r);
}

Eigen::Matrix3cd green_dyadic(const Eigen::Vector3d &rv, double k) {
  double r = rv.norm();
  Eigen::Vector3d rhat = rv / r;

  complex Phi = green_scalar(r, k);

  Eigen::Matrix3d Id = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d rr = rhat * rhat.transpose();

  complex traverse = k * k + 1i * k / r - 1 / (r * r);
  complex radial = -k * k - 3i * k / r + 3 / (r * r);

  return (1i / k) * Phi * (traverse * Id + radial * rr);
}

Eigen::Vector3cd incident_field(
  const Eigen::Vector3d &x,
  const Eigen::Vector3d &z,
  double k,
  const Eigen::Vector3d &p
) {
  return green_dyadic(x - z, k) * p;
}

// TODO: is this a good choice?
// Fibonacci sphere grid (golden-ratio spiral) on radius R
Eigen::MatrixXd sphere_grid(double R, int N) {
  double golden = (1.0 + std::sqrt(5.0)) / 2.0;
  Eigen::MatrixXd pts(N, 3);
  for (int i = 0; i < N; ++i) {
    double phi = std::acos(1.0 - 2.0 * (i + 0.5) / N);
    double theta = 2.0 * PI * (i + 0.5) / golden;
    pts.row(i) = Eigen::Vector3d(R * std::cos(theta) * std::sin(phi),
                                 R * std::sin(theta) * std::sin(phi),
                                 R * std::cos(phi));
  }
  return pts;
}

// Build a tangential basis {e1, e2} at a point on the unit sphere
// {e1,e2,nhat} form a right-handed system
void tangential_basis(
  const Eigen::Vector3d &nhat,
  Eigen::Vector3d &e1,
  Eigen::Vector3d &e2
) {
  Eigen::Vector3d ref = (std::abs(nhat(0)) < 0.9) ?
    Eigen::Vector3d(1, 0, 0) : Eigen::Vector3d(0, 1, 0);
  e1 = ref.cross(nhat).normalized(); // a vector perpendicular to nhat
  e2 = nhat.cross(e1).normalized(); // a vector perpendicular to both nhat and e1
}

// Scale a sphere geometry to an ellipsoid by modifying control points.
// The NURBS data is stored in projective coordinates (wx, wy, wz, w),
// so scaling physical coords by (a, b, c) means multiplying the first
// three components by (a, b, c) while leaving w unchanged.
Bembel::Geometry make_ellipsoid(
  const std::string &sphere_dat_path, double a,
  double b, double c
) {
  auto patches = Bembel::LoadGeometryFile(sphere_dat_path);
  for (auto &patch : patches) {
    int n = patch.data_.size() / 4;
    for (int i = 0; i < n; ++i) {
      patch.data_[4 * i + 0] *= a;
      patch.data_[4 * i + 1] *= b;
      patch.data_[4 * i + 2] *= c;
    }
  }
  return Bembel::Geometry(patches);
}

// Write helpers
void write_complex_matrix(const std::string &path, const Eigen::MatrixXcd &mat) {
  std::ofstream f(path);
  f << std::setprecision(15);
  for (int i = 0; i < mat.rows(); ++i) {
    for (int j = 0; j < mat.cols(); ++j) {
      if (j > 0) f << " ";
      f << mat(i, j).real() << " " << mat(i, j).imag();
    }
    f << "\n";
  }
}

void write_real_matrix(const std::string &path, const Eigen::MatrixXd &mat) {
  std::ofstream f(path);
  f << std::setprecision(15);
  for (int i = 0; i < mat.rows(); ++i) {
    for (int j = 0; j < mat.cols(); ++j) {
      if (j > 0) f << " ";
      f << mat(i, j);
    }
    f << "\n";
  }
}

int main(int argc, char *argv[]) {
  using namespace Bembel;
  using namespace Eigen;

  // --- parameters ---
  std::string sphere_dat = "res/sphere.dat";
  double k = 2.0;
  int refinement = 1;//3;
  int poly_deg = 1;//2;
  int n_tx = 50;
  int n_rx = 200;

  if (argc > 1) sphere_dat = argv[1];
  if (argc > 2) k = std::stod(argv[2]);
  if (argc > 3) refinement = std::stoi(argv[3]);
  if (argc > 4) poly_deg = std::stoi(argv[4]);
  if (argc > 5) n_tx = std::stoi(argv[5]);
  if (argc > 6) n_rx = std::stoi(argv[6]);

  // Bembel's kernel is exp(-i * wavenumber * r) / (4 pi r).
  // Setting wavenumber = -k gives exp(+ikr) / (4 pi r),
  // matching the e^{+ikr} convention used in the incident field.
  complex wavenumber(-k, 0.0);

  std::cout << "=== PEC Ellipsoidal Cavity BEM Solver ===" << std::endl;
  std::cout << "  sphere.dat:    " << sphere_dat << std::endl;
  std::cout << "  k:             " << k << std::endl;
  std::cout << "  refinement:    " << refinement << std::endl;
  std::cout << "  poly degree:   " << poly_deg << std::endl;
  std::cout << "  transmitters:  " << n_tx << std::endl;
  std::cout << "  receivers:     " << n_rx << std::endl;

  // --- geometry: ellipsoid with semi-axes (4, 4, 6) ---
  Geometry geometry = make_ellipsoid(sphere_dat, 4.0, 4.0, 6.0);

  // --- transmitter and receiver grids on Lambda (unit sphere) ---
  MatrixXd tx_pts = sphere_grid(1.0, n_tx);
  MatrixXd rx_pts = sphere_grid(1.0, n_rx);

  // build dipole configs: 2 tangential polarizations per transmitter
  struct Dipole {
    Vector3d z, p;
  };
  std::vector<Dipole> dipoles;
  for (int i = 0; i < n_tx; ++i) {
    Vector3d nhat = tx_pts.row(i).normalized();
    Vector3d e1, e2;
    tangential_basis(nhat, e1, e2);
    dipoles.push_back({tx_pts.row(i), e1});
    dipoles.push_back({tx_pts.row(i), e2});
  }
  int n_dipoles = dipoles.size();

  std::cout << "  dipoles:       " << n_dipoles << std::endl;

  // --- assemble BEM operator (done once) ---
  std::cout << "\nAssembling operator..." << std::flush;
  AnsatzSpace<MaxwellSingleLayerOperator> ansatz_space(geometry, refinement, poly_deg);

  DiscreteOperator<H2Matrix<complex>, MaxwellSingleLayerOperator> disc_op(ansatz_space);
  disc_op.get_linear_operator().set_wavenumber(wavenumber);
  disc_op.compute();
  std::cout << " done." << std::endl;

  GMRES<H2Matrix<complex>, IdentityPreconditioner> gmres;
  gmres.compute(disc_op.get_discrete_operator());
  gmres.set_restart(2000);
  gmres.setTolerance(1e-8);

  DiscretePotential<MaxwellSingleLayerPotential<MaxwellSingleLayerOperator>,
                    MaxwellSingleLayerOperator>
      disc_pot(ansatz_space);
  disc_pot.get_potential().set_wavenumber(wavenumber);

  // --- solve for each dipole excitation ---
  // Es_all(r, 3*d + c) = c-th component of E^s at receiver r for dipole d
  MatrixXcd Es_all(n_rx, 3 * n_dipoles);

  for (int d = 0; d < n_dipoles; ++d) {
    Vector3d z = dipoles[d].z;
    Vector3d p = dipoles[d].p;

    std::cout << "  dipole " << d + 1 << "/" << n_dipoles << std::flush;

    // RHS: pass -E^i to RotatedTangentialTrace
    const std::function<VectorXcd(Vector3d)> neg_Ei =
        [k, z, p](Vector3d x) -> VectorXcd {
      return -incident_field(x, z, k, p);
    };

    DiscreteLinearForm<RotatedTangentialTrace<complex>,
                       MaxwellSingleLayerOperator>
        disc_lf(ansatz_space);
    disc_lf.get_linear_form().set_function(neg_Ei);
    disc_lf.compute();

    VectorXcd rho = gmres.solve(disc_lf.get_discrete_linear_form());
    std::cout << "  GMRES iters: " << gmres.iterations()
              << "  residual: " << gmres.error() << std::endl;

    disc_pot.set_cauchy_data(rho);
    
    MatrixXcd Es = disc_pot.evaluate(rx_pts);  // (n_rx x 3)

    Es_all.block(0, 3 * d, n_rx, 3) = Es;

    double max_norm = Es.rowwise().norm().maxCoeff();
    std::cout << "  |Es|_max = " << max_norm << std::endl;
  }

  // --- assemble tangential response operator T ---
  // T(2*r + c, d) = projection of E^s_d(x_r) onto c-th tangential basis
  MatrixXcd T(2 * n_rx, n_dipoles);
  for (int r = 0; r < n_rx; ++r) {
    Vector3d nhat = rx_pts.row(r).normalized();
    Vector3d e1, e2;
    tangential_basis(nhat, e1, e2);
    Vector3cd e1c = e1.cast<complex>();
    Vector3cd e2c = e2.cast<complex>();

    for (int d = 0; d < n_dipoles; ++d) {
      Vector3cd Es_rd = Es_all.block(r, 3 * d, 1, 3).transpose();
      T(2 * r, d) = Es_rd.dot(e1c);
      T(2 * r + 1, d) = Es_rd.dot(e2c);
    }
  }

  // --- write output ---
  write_real_matrix("out/receivers.dat", rx_pts);
  write_complex_matrix("out/T_matrix.dat", T);
  write_complex_matrix("out/Es_all.dat", Es_all);

  {
    std::ofstream f("out/dipoles.dat");
    f << std::setprecision(15);
    for (auto &dip : dipoles)
      f << dip.z.transpose() << " " << dip.p.transpose() << "\n";
  }

  std::cout << "\n=== Done ===" << std::endl;
  std::cout << "Output:" << std::endl;
  std::cout << "  receivers.dat  (" << n_rx << " x 3)" << std::endl;
  std::cout << "  dipoles.dat    (" << n_dipoles << " x 6)" << std::endl;
  std::cout << "  Es_all.dat     (" << n_rx << " x " << 3 * n_dipoles
            << ") [re im ...]" << std::endl;
  std::cout << "  T_matrix.dat   (" << 2 * n_rx << " x " << n_dipoles
            << ") [re im ...]" << std::endl;


  return 0;
}


