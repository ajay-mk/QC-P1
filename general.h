#include <iostream>
#include <string>
#include <Eigen/Eigenvalues>

int count_lines(std::string filename);
auto Diagonalize(Eigen::MatrixXd M);
void Diag2(Eigen::MatrixXd M, Eigen::MatrixXd evals, Eigen::MatrixXd evecs);
