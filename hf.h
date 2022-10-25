#include <iostream>
#include <string>
#include "Eigen/Eigen"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"


// Functions for reading input
double read_enuc(std::string int_path);
Eigen::MatrixXd read_1e_ints(std::string int_path, std::string int_file);
Eigen::MatrixXd read_2e_ints(std::string int_path, std::string int_file);

// Functions
struct diag_results {Eigen::MatrixXd evecs; Eigen::MatrixXd evals;};
diag_results Diag_M(Eigen::MatrixXd M);
diag_results Diag(Eigen::MatrixXd M);
Eigen::MatrixXd build_density(Eigen::MatrixXd Coeff, int nocc);
double scf_energy(const Eigen::MatrixXd& Density, const Eigen::MatrixXd& H, const Eigen::MatrixXd& Fock);
Eigen::MatrixXd fock_build(Eigen::MatrixXd H, Eigen::MatrixXd Density, Eigen::MatrixXd two_e);