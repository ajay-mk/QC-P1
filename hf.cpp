// Hartree-Fock algorithm for practice
// Integrals are read from external files in this version
#include "hf.h"
#include "molecule.h"
#include <string>
#include <fstream>
#include "Eigen/Eigen"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

using std::cout, std::endl;

// Since we are considering water molecule in STO-3G basis
const int nao = 7;

// Compound Indices
int CompoundIndex(int a, int b) {
    int ab;
    if (a > b)
        ab = (a * (a + 1) / 2) + b;
    else
        ab = (b * (b + 1) / 2) + a;
    return ab;
}
// Main Hartree-Fock Function
hartree_fock::hartree_fock(std::string filename, std::string int_path){
    Molecule mol(filename);
    mol.print_info();
    mol.print_geometry();

    // Reading Nuclear Repulsion
    double enuc  = read_enuc(int_path);
    cout << "----Nuclear repulsion energy----" << endl << enuc  << endl << endl;
    // Reading 1e integrals
    Eigen::Matrix<double, nao, nao> S = read_1e_ints(int_path, "/overlap.dat");
    cout << "----Overlap Integrals----" << endl << S << endl;
    cout << endl;
    Eigen::Matrix<double, nao, nao> T = read_1e_ints(int_path, "/ke.dat");
    cout << "----Kinetic Energy Integrals----" << endl << T << endl;
    cout << endl;
    Eigen::Matrix<double, nao, nao> V = read_1e_ints(int_path, "/v.dat");
    cout << "----Nuclear Attraction Integrals----" << endl << V << endl;
    cout << endl;
    // Forming Core Hamiltonian
    Eigen::Matrix<double, nao, nao> H = T + V;
    cout << "----Core Hamiltonian----" << endl << H << endl;
    cout << endl;

    // Reading 2-e integral
    Eigen::MatrixXd two_e = read_2e_ints(int_path, "/eri.dat");

    // Diagonalizing Overlap Matrix
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver1(S);
    auto S_evals_c = solver1.eigenvalues();
    auto S_evecs = solver1.eigenvectors(); //lambda_S

    // S_evals is a column matrix, we need to build a diagonal matrix with its elements
    Eigen::MatrixXd S_evals = Eigen::MatrixXd ::Zero(S.rows(), S.rows()); //L_S
    for (int i =0; i < S.rows(); i++)
    {
        for (int j = 0; j < S.rows(); j++)
        {
            if(i==j)
                S_evals(i, i) = S_evals_c(i);
        }
    }

    // Building Symmetric Orthogonolization Matrix
    for (int i=0; i < S.rows(); i++)
    {
        S_evals(i, i) = pow(S_evals(i, i), -0.5); // Square root of S_evals
    }
    auto SOM = S_evecs * S_evals * S_evecs.transpose();
    cout << "----Symmetric Orthogonalization Matrix----" << endl;
    cout << SOM << endl;

    cout << endl;

    // Building Initial Fock Matrix - Using Core Hamiltonian
    auto F0 = SOM.transpose() * H * SOM ;
    cout << "----Initial Fock Matrix----" << endl;
    cout << F0 << endl;

    // Diagonalizing Initial Fock Matrix
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver2(F0);
    auto C0_ = solver2.eigenvectors();
    auto eps = solver2.eigenvalues(); // Initial Orbital Energy

    // Transformation of Eigenvectors into original AO basis
    auto C0 = SOM * C0_;
    cout << endl;
    cout << "----Initial Coefficient Matrix----" << endl;
    cout << C0 << endl;


}

double hartree_fock::read_enuc(std::string int_path)
// Reading nuclear repulsion energy
// Expecting 'enuc.dat' at int_path
{
    std::string file_path = int_path + "/enuc.dat";
    std::ifstream enuc_file(file_path);
    double enuc;
    if (enuc_file.is_open())
        enuc_file >> enuc;
    else
        cout << "Unable to open enuc.dat" << endl;
    return enuc;
}

Eigen::MatrixXd hartree_fock::read_1e_ints(std::string int_path, std::string int_file)
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
Eigen::MatrixXd hartree_fock::read_2e_ints(std::string int_path, std::string int_file)
{
    std::string filename = int_path + int_file;
    // Here we use the compound indices defined earlier
    // What will be the size of the matrix? num = M(M+1)/2; M = nao(nao+1)/2
    int M = nao * (nao + 1)/2;
    int num = M * (M +1)/2;
    Eigen::MatrixXd two_e = Eigen::MatrixXd::Zero(num, 1);
    std::ifstream int_input(filename);
    while(!int_input.eof())
    {
        int i, j, k, l, ij, kl;
        double temp;
        int_input >> i >> j >> k >> l >> temp;
        ij = CompoundIndex(i-1, j-1);
        kl = CompoundIndex(k-1, l-1);
        two_e (CompoundIndex(ij, kl), 0) = temp;
    }
    return two_e;
}
