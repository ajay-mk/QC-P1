// Hartree-Fock algorithm for practice
// Integrals are read from external files in this version
#include "hf.h"
#include "molecule.h"
#include <string>
#include <fstream>
#include <array>
#include "Eigen/Eigen"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

using std::cout, std::endl;

// Since we are considering water molecule in STO-3G basis
const int nao = 7;

// Compound Indices
int CompoundIndex(int i, int j) {
    int ij;
    if (i > j)
        ij = i * (i + 1) / (2 + j);
    else
        ij = j * (j + 1) / (2 + i);
    return ij;
}

hartree_fock::hartree_fock(std::string filename, std::string int_path){
    Molecule mol(filename);
    mol.print_info();
    mol.print_geometry();

    // Reading Nuclear Repulsion
    double enuc  = read_enuc(int_path);
    cout << "Nuclear repulsion energy: " << enuc  << endl << endl;
    // Reading 1e integrals
    Eigen::Matrix<double, nao, nao> S = read_1e_ints(int_path, "/overlap.dat");
    cout << "Overlap Integrals" << endl << S << endl;
    cout << endl;
    Eigen::Matrix<double, nao, nao> T = read_1e_ints(int_path, "/ke.dat");
    cout << "Kinetic Energy Integrals" << endl << T << endl;
    cout << endl;
    Eigen::Matrix<double, nao, nao> V = read_1e_ints(int_path, "/v.dat");
    cout << "Nuclear Attraction Integral" << endl << V << endl;
    cout << endl;
    // Forming Core Hamiltonian
    Eigen::Matrix<double, nao, nao> H = T + V;
    cout << "Core Hamiltonian" << endl << H << endl;
    cout << endl;

    // Reading 2-e integral
    Eigen::MatrixXd two_e = read_2e_ints(int_path, "/eri.dat");
    cout << two_e << endl;
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
    int num = (((nao * (nao + 1) / 2) + 1) * (nao * (nao + 1) / 2)/2);
    Eigen::MatrixXd two_e = Eigen::MatrixXd::Zero(1, num);
    std::ifstream int_input(filename);
    while(!int_input.eof())
    {
        int mu, nu, lambda, sigma, mu_nu, lambda_sigma;
        double temp;
        int_input >> mu >> nu >> lambda >> sigma >> temp;
        mu_nu = CompoundIndex(mu, nu);
        lambda_sigma = CompoundIndex(lambda, sigma);
        int mu_nu_lambda_sigma = CompoundIndex(mu_nu, lambda_sigma);
        two_e (0 , mu_nu_lambda_sigma) = temp;
    }
    return two_e;
}
