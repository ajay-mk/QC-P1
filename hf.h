#include <iostream>
#include <string>
#include "Eigen/Eigen"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"


// Functions for reading input
double read_enuc(const std::string& int_path);
Eigen::MatrixXd read_1e_ints(const std::string& int_path, const std::string& int_file);
Eigen::MatrixXd read_2e_ints(const std::string& int_path, const std::string& int_file);

// Functions
struct diag_results {Eigen::MatrixXd evecs; Eigen::MatrixXd evals;};
diag_results Diag_M(const Eigen::MatrixXd& M);
diag_results Diag(const Eigen::MatrixXd& M);
Eigen::MatrixXd build_density(const Eigen::MatrixXd& Coeff, int nocc);
double scf_energy(const Eigen::MatrixXd& Density, const Eigen::MatrixXd& H, const Eigen::MatrixXd& Fock);
Eigen::MatrixXd fock_build(const Eigen::MatrixXd& H, const Eigen::MatrixXd& Density, const Eigen::MatrixXd& two_e);