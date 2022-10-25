#include "hf.h"
#include "molecule.h"
#include <string>
#include <fstream>
#include "Eigen/Eigenvalues"

using std::cout, std::endl;

int nao = 7;

// Compound Indices
int CompoundIndex(int a, int b, int c, int d) {
    double ab;
    if (a > b)
        ab = a * (a + 1)/2 + b;
    else
        ab = b * (b + 1)/2 + a;
    double cd;
    if (c > d)
        cd = c * (c + 1)/2 + d;
    else
        cd = d * (d + 1)/2 + c;

    double abcd;
    if (ab > cd)
        abcd = ab * (ab + 1)/2 + cd;
    else
        abcd = cd * (cd + 1)/2 + ab;
    return int(abcd);

}

// Function Definitions

double read_enuc(std::string int_path)
// Reading nuclear repulsion energy
// Expecting 'enuc.dat' at int_path
{
    std::string file_path = int_path + "/enuc.dat";
    std::ifstream enuc_file(file_path);
    double num;
    if (enuc_file.is_open())
        enuc_file >> num;
    else
        cout << "Unable to open enuc.dat" << endl;
    return num;
}

Eigen::MatrixXd read_1e_ints(std::string int_path, std::string int_file)
{
    Eigen::MatrixXd one_e_matrix = Eigen::MatrixXd::Zero(nao, nao);
    std::string filename = int_path + int_file;
    std::ifstream int_input(filename);
    while(!int_input.eof()) {
        int a, b; double temp;
        int_input >> a >> b >> temp;
        one_e_matrix(a-1, b-1) = temp;
        one_e_matrix(b-1, a-1) = temp;
    }
    return one_e_matrix;

}
// We have four-dimensional matrix of 2e-integrals, because of the eight-fold symmetry
// we only need the lower triangle and this can be stored into a column matrix
Eigen::MatrixXd read_2e_ints(std::string int_path, std::string int_file)
{
    std::string filename = int_path + int_file;
    // Here we use the compound indices defined earlier
    // What will be the size of the matrix? num = M(M+1)/2; M = nao(nao+1)/2
//    int M = nao * (nao + 1)/2;
//    int num = M * (M +1)/2;
    int num = CompoundIndex(nao, nao, nao, nao);
    Eigen::MatrixXd two_e = Eigen::MatrixXd::Zero(num, 1);
    std::ifstream int_input(filename);
    while(!int_input.eof())
    {
        int i, j, k, l;
        double temp;
        int_input >> i >> j >> k >> l >> temp;
        two_e (CompoundIndex(i-1, j-1, k-1, l-1), 0) = temp;
    }
    return two_e;
}
diag_results Diag(Eigen::MatrixXd M) {
    diag_results results;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(M);
    results.evecs = solver.eigenvectors();
    results.evals = solver.eigenvalues();
    return results;
}


diag_results Diag_M(Eigen::MatrixXd M) {
    int num = M.rows();
    diag_results results;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(M);
    results.evecs = solver.eigenvectors();
    auto evals_c= solver.eigenvalues();
    results.evals = Eigen::MatrixXd::Zero(num, num);
    for(int i=0; i < num; i++){
        for(int j=0; j < num; j++){
            if(i==j)
                results.evals(i, i) = evals_c(i);
            else
                results.evals(i, j) = 0;
        }
    }
    return results;
}

Eigen::MatrixXd build_density(Eigen::MatrixXd Coeff, int nocc) {
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(nao, nao);
    for(int i = 0; i < nao; i++){
        for(int j =0; j < nao; j++){
            // Sum over occupied orbitals
            for(int m =0; m < nocc; m++){
                D(i, j) += Coeff(i, m) * Coeff(j, m);
            }
        }
    }
    return D;
}

double scf_energy(const Eigen::MatrixXd& Density, const Eigen::MatrixXd& H, const Eigen::MatrixXd& Fock)
{
    double e_hf = 0;
    for(int i=0; i < nao; i++){
        for(int j=0; j < nao; j++){
            e_hf += Density(i, j) * (H(i, j) + Fock (i, j));
        }
    }
    return e_hf;
}

Eigen::MatrixXd fock_build(Eigen::MatrixXd H, Eigen::MatrixXd Density, Eigen::MatrixXd two_e)
{
    auto F = H;
    for (int i = 0; i < nao; i++) {
        for (int j = 0; j < nao; j++) {
            for (int k = 0; k < nao; k++) {
                for (int l = 0; l < nao; l++) {
                    F(i, j) += Density(k, l) * (2 * two_e(CompoundIndex(i, j, k, l), 0) - two_e(CompoundIndex(i, k, j, l), 0));
                }
            }
        }
    }
    return F;
}