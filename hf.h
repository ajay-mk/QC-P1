#include <iostream>
#include <string>
#include "Eigen/Eigen"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"


class hartree_fock{
public:
    hartree_fock(std::string filename, std::string int_path);
    double read_enuc(std::string int_path);
    Eigen::MatrixXd read_1e_ints(std::string int_path, std::string int_file);
    Eigen::MatrixXd read_2e_ints(std::string int_path, std::string int_file);
    Eigen::MatrixXd S;


private:
    double nelectrons;

};
