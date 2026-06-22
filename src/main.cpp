#include <Bembel/AnsatzSpace>
#include <Bembel/Geometry>
#include <Bembel/IO>
#include <Bembel/LinearForm>
#include <Bembel/Maxwell>
#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <filesystem>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

static constexpr double PI = 3.14159265358979323846;

using namespace std::complex_literals;
using complex = std::complex<double>;

struct Config {
  double k;
  Eigen::Vector3d semiaxes;
  Eigen::MatrixXd points;
  Eigen::MatrixXd e1;
  Eigen::MatrixXd e2;
};

Config load_config(const std::string &path) {
  std::ifstream f(path);
  std::string line;
  Config s;
  std::vector<std::array<double, 9>> rows;
  bool header_read = false;

  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream iss(line);
    if (!header_read) {
      double a, b, c;
      int n;
      iss >> s.k >> a >> b >> c >> n;
      s.semiaxes = Eigen::Vector3d(a, b, c);
      header_read = true;
    } else {
      std::array<double, 9> r;
      for (double &v : r) iss >> v;
      rows.push_back(r);
    }
  }

  int n = rows.size();
  s.points.resize(n, 3);
  s.e1.resize(n, 3);
  s.e2.resize(n, 3);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < 3; ++j) {
      s.points(i, j) = rows[i][j];
      s.e1(i, j) = rows[i][3 + j];
      s.e2(i, j) = rows[i][6 + j];
    }
  return s;
}

complex green_scalar(double r, double k) {
  return std::exp(1i * k * r) / (4 * PI * r);
}

Eigen::Matrix3cd green_dyadic(const Eigen::Vector3d &rv, double k) {
  double r = rv.norm();
  Eigen::Vector3d rhat = rv / r;
  complex Phi = green_scalar(r, k);
  Eigen::Matrix3d Id = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d rr = rhat * rhat.transpose();
  complex transverse = k * k + 1i * k / r - 1 / (r * r);
  complex radial = -k * k - 3i * k / r + 3 / (r * r);
  return (1i / k) * Phi * (transverse * Id + radial * rr);
}

Eigen::Vector3cd incident_field(const Eigen::Vector3d &x, const Eigen::Vector3d &z,
                                double k, const Eigen::Vector3d &p) {
  return green_dyadic(x - z, k) * p;
}

Bembel::Geometry make_ellipsoid(const std::string &path, double a, double b, double c) {
  auto patches = Bembel::LoadGeometryFile(path);
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

void write_complex_matrix(const std::string &path, const Eigen::MatrixXcd &mat) {
  std::ofstream f(path);
  if (!f) throw std::runtime_error("cannot open " + path);
  f << std::setprecision(15);
  for (int i = 0; i < mat.rows(); ++i) {
    for (int j = 0; j < mat.cols(); ++j) {
      if (j > 0) f << " ";
      f << mat(i, j).real() << " " << mat(i, j).imag();
    }
    f << "\n";
  }
}

// ksweep: assemble the EFIE operator over a wavenumber grid and record its
// condition number per k. No RHS, solve, or operator T -- the conditioning of
// the boundary integral operator alone detects interior cavity resonances.
int run_ksweep(const Config &config, const std::string &sphere_dat,
               int refinement, int poly_deg, double kmin, double kmax, int nk) {
  using namespace Bembel;
  using namespace Eigen;

  Geometry geometry = make_ellipsoid(sphere_dat, config.semiaxes(0),
                                     config.semiaxes(1), config.semiaxes(2));
  AnsatzSpace<MaxwellSingleLayerOperator> ansatz_space(geometry, refinement, poly_deg);

  std::cout << "=== BEM wavenumber sweep ===\n";
  std::cout << "  refinement:  " << refinement << "\n";
  std::cout << "  poly degree: " << poly_deg << "\n";
  std::cout << "  dofs:        " << ansatz_space.get_number_of_dofs() << "\n";
  std::cout << "  k range:     [" << kmin << ", " << kmax << "], " << nk << " points\n"
            << std::flush;

  std::ofstream out("out/ksweep.csv");
  if (!out) throw std::runtime_error("cannot open out/ksweep.csv");
  out << "k,cond\n" << std::setprecision(10);
  for (int i = 0; i < nk; ++i) {
    double kk = nk > 1 ? kmin + (kmax - kmin) * i / (nk - 1) : kmin;
    DiscreteOperator<MatrixXcd, MaxwellSingleLayerOperator> disc_op(ansatz_space);
    disc_op.get_linear_operator().set_wavenumber(complex(kk, 0.0));
    disc_op.compute();
    double cond = 1.0 / PartialPivLU<MatrixXcd>(disc_op.get_discrete_operator()).rcond();
    out << kk << "," << cond << "\n" << std::flush;
    std::cout << "k=" << kk << " cond=" << cond << "\n" << std::flush;
  }
  std::cout << "wrote out/ksweep.csv\n";
  return 0;
}

// operator: assemble the dipole reaction operator T at a fixed (refinement,
// poly_deg) and write it to disk, with the EFIE system condition number.
int run_operator(const Config &config, const std::string &sphere_dat,
                 int refinement, int poly_deg) {
  using namespace Bembel;
  using namespace Eigen;

  double k = config.k;
  complex wavenumber(k, 0.0);
  int n_points = config.points.rows();
  int n_dipoles = 2 * n_points;

  Geometry geometry = make_ellipsoid(sphere_dat, config.semiaxes(0),
                                     config.semiaxes(1), config.semiaxes(2));

  std::vector<Vector3d> z(n_dipoles), pol(n_dipoles);
  for (int i = 0; i < n_points; ++i) {
    z[2 * i] = z[2 * i + 1] = config.points.row(i);
    pol[2 * i] = config.e1.row(i);
    pol[2 * i + 1] = config.e2.row(i);
  }

  AnsatzSpace<MaxwellSingleLayerOperator> ansatz_space(geometry, refinement, poly_deg);
  const int N = ansatz_space.get_number_of_dofs();

  std::cout << "=== PEC Ellipsoidal Cavity BEM Solver ===\n";
  std::cout << "  k:           " << k << "\n";
  std::cout << "  refinement:  " << refinement << "\n";
  std::cout << "  poly degree: " << poly_deg << "\n";
  std::cout << "  dipoles:     " << n_dipoles << "\n";
  std::cout << "  dofs:        " << N << "\n";
  std::cout << "  dense matrix: ~" << 16.0 * double(N) * double(N) / 1e9 << " GB\n"
            << std::flush;

  MatrixXcd B(N, n_dipoles);
  for (int d = 0; d < n_dipoles; ++d) {
    const std::function<VectorXcd(Vector3d)> neg_Ei =
        [&, d](Vector3d x) -> VectorXcd { return -incident_field(x, z[d], k, pol[d]); };
    DiscreteLinearForm<RotatedTangentialTrace<complex>, MaxwellSingleLayerOperator>
        disc_lf(ansatz_space);
    disc_lf.get_linear_form().set_function(neg_Ei);
    disc_lf.compute();
    B.col(d) = disc_lf.get_discrete_linear_form();
  }

  std::cout << "\nAssembling operator (dense) and solving..." << std::flush;

  // The Maxwell single-layer (EFIE) is first-kind and badly conditioned, so an
  // unpreconditioned iterative solve is hopeless. Assemble the dense operator,
  // factor once, and direct-solve all dipole RHS at once.
  DiscreteOperator<MatrixXcd, MaxwellSingleLayerOperator> disc_op(ansatz_space);
  disc_op.get_linear_operator().set_wavenumber(wavenumber);
  disc_op.compute();
  PartialPivLU<MatrixXcd> lu(disc_op.get_discrete_operator());
  // Cheap 1-norm reciprocal condition estimate of the EFIE system matrix (reuses
  // the LU; no extra factorization), reported alongside EP-GP's system cond.
  double cond_sys = 1.0 / lu.rcond();
  std::cout << "cond(A) = " << cond_sys << "\n" << std::flush;
  MatrixXcd Rho = lu.solve(B);

  DiscretePotential<MaxwellSingleLayerPotential<MaxwellSingleLayerOperator>,
                    MaxwellSingleLayerOperator>
      disc_pot(ansatz_space);
  disc_pot.get_potential().set_wavenumber(wavenumber);

  MatrixXcd Es_all(n_points, 3 * n_dipoles);

  for (int d = 0; d < n_dipoles; ++d) {
    disc_pot.set_cauchy_data(Rho.col(d));
    Es_all.block(0, 3 * d, n_points, 3) = disc_pot.evaluate(config.points);
  }


  // Bembel works in the e^{-iωt} convention; conjugate to match e^{+iωt}.
  Es_all = Es_all.conjugate();

  write_complex_matrix("out/Es_all.dat", Es_all);

  MatrixXcd T(n_dipoles, n_dipoles);
  for (int r = 0; r < n_points; ++r) {
    Vector3cd nhat = config.points.row(r).normalized().cast<complex>();
    Vector3cd e1c = config.e1.row(r).cast<complex>();
    Vector3cd e2c = config.e2.row(r).cast<complex>();
    for (int d = 0; d < n_dipoles; ++d) {
      Vector3cd Es = Es_all.block(r, 3 * d, 1, 3).transpose();
      Vector3cd Es_t = Es - Es.dot(nhat) * nhat;
      T(2 * r, d) = Es_t.dot(e1c);
      T(2 * r + 1, d) = Es_t.dot(e2c);
    }
  }

  double norm_T = T.norm();
  double rel_asym = (T - T.transpose()).norm() / norm_T;
  std::cout << "\n||T||_F = " << norm_T << "\n";
  std::cout << "||T - T^T||_F / ||T||_F = " << rel_asym << "\n";

  write_complex_matrix("out/T_matrix.dat", T);

  std::cout << "wrote out/T_matrix.dat (" << n_dipoles << " x " << n_dipoles << ")\n";

  return 0;
}

int main(int argc, char *argv[]) {
  // usage:
  //   my_project operator <config> <refinement> <poly_deg>
  //   my_project ksweep   <config> <refinement> <poly_deg> <kmin> <kmax> <nk>
  if (argc < 2) {
    std::cerr << "usage: " << argv[0]
              << " operator|ksweep <config> <refinement> <poly_deg> [kmin kmax nk]\n";
    return 1;
  }
  std::string mode = argv[1];
  std::string config_path = argc > 2 ? argv[2] : "res/config_ellipse.txt";
  int refinement = argc > 3 ? std::stoi(argv[3]) : 1;
  int poly_deg = argc > 4 ? std::stoi(argv[4]) : 2;

  std::filesystem::create_directories("out");
  Config config = load_config(config_path);
  // sphere.dat lives next to the config, so it resolves regardless of cwd.
  std::string sphere_dat =
      (std::filesystem::path(config_path).parent_path() / "sphere.dat").string();

  if (mode == "operator") {
    return run_operator(config, sphere_dat, refinement, poly_deg);
  } else if (mode == "ksweep") {
    if (argc < 8) {
      std::cerr << "ksweep needs: <config> <refinement> <poly_deg> <kmin> <kmax> <nk>\n";
      return 1;
    }
    return run_ksweep(config, sphere_dat, refinement, poly_deg,
                      std::stod(argv[5]), std::stod(argv[6]), std::stoi(argv[7]));
  }
  std::cerr << "unknown mode '" << mode << "' (expected operator|ksweep)\n";
  return 1;
}
