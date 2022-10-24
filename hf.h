#include <iostream>
#include <string>
#include "general.h"
#include "Eigen/Eigen"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"


class hartree_fock{
public:
    hartree_fock(std::string filename, std::string int_path);
    double read_enuc(std::string int_path);
    // Functions for reading input
    Eigen::MatrixXd read_1e_ints(std::string int_path, std::string int_file);
    Eigen::MatrixXd read_2e_ints(std::string int_path, std::string int_file);
    double enuc;
    // Matrices for integrals
    Eigen::MatrixXd S;
    Eigen::MatrixXd T;
    Eigen::MatrixXd V;
    Eigen::MatrixXd two_e;
    // Matrices for SCF loop
    Eigen::MatrixXd H;
    double hf_energy;

    // Functions
    struct diag_results{Eigen::MatrixXd evecs; Eigen::MatrixXd evals;};
    diag_results Diag_M(Eigen::MatrixXd M, diag_results);
    diag_results Diag(Eigen::MatrixXd M, diag_results);
    Eigen::MatrixXd build_density(Eigen::MatrixXd Coeff, int nocc);
    double scf_energy(Eigen::MatrixXd D, Eigen::MatrixXd H, Eigen::MatrixXd F);
    Eigen::MatrixXd fock_build(Eigen::MatrixXd H, Eigen::MatrixXd D, Eigen::MatrixXd two_e);

private:
    double nelectrons;

};
