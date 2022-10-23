#include "general.h"
#include <fstream>

using std::cout;
using std::endl;

auto Diagonalize(Eigen::MatrixXd M){
    int num_rows = M.rows();
    Eigen::MatrixXd M_evals = Eigen::MatrixXd::Zero(M.rows(), M.rows());
    struct result{Eigen::MatrixXd evecs, evals;};
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(M);
    auto M_evecs = solver.eigenvectors();
    auto M_evals_vec = solver.eigenvalues();
    for(int i=0; i < num_rows; i++){
        for (int j=0; j< num_rows; j++){
            if(i==j)
                M_evals(i, i) = M_evals_vec(i);
        }
    }
    return result {M_evecs, M_evals};
}

void Diag2(Eigen::MatrixXd M, Eigen::MatrixXd evals, Eigen::MatrixXd evecs) {
    const int num_rows = M.rows();
    evecs.resize(num_rows, num_rows);
    evals.resize(num_rows, num_rows);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(M);
    evecs = solver.eigenvectors();
    auto evals_vec = solver.eigenvalues();
    for (int i=0; i<num_rows; i++) {
        for (int j=0; j<num_rows; j++) {
            if (i==j) {
                (evals)(i,i) = evals_vec(i);
            } else {
                (evals)(i,j) = 0;
            }
        }
    }
}

// Function for counting number of lines in a file
int count_lines(std::string filename)
{
    int nlines;
    std::string line;
    std::ifstream input (filename);
    while(!input.eof())
        ++ nlines;
    return nlines;
}