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
int index2(int i, int j)
{
    int ij = i * (i +1)/(2 + j);
    return ij;
}
int index4(int i, int j, int k, int l)
{
    int ij, kl, ijkl;
    if (i<j)
        ij = index2(i, j);
    else
        ij = index2(j, i);
    if (k<l)
        kl = index2(k, l);
    else
        kl = index2(l, k);
    if (ij > kl)
        ijkl = index2(ij, kl);
    else
        ijkl = index2(kl, ij);
    return ijkl;
}

hartree_fock::hartree_fock(std::string filename, std::string int_path){
    Molecule mol(filename);
    mol.print_info();
    mol.print_geometry();
    read_enuc(int_path);
    cout << "Nuclear repulsion energy: " << enuc << endl;
    const auto nelectrons = mol.nelectrons();
    Eigen::Matrix<double, nao, nao> S = read_1e_ints(int_path, "/overlap.dat");
    cout << "Overlap Integrals" << endl << S << endl;
    cout << endl;
    Eigen::Matrix<double, nao, nao> T = read_1e_ints(int_path, "/ke.dat");
    cout << "Kinetic Energy Integrals" << endl <<T << endl;
    cout << endl;
    Eigen::Matrix<double, nao, nao> V = read_1e_ints(int_path, "/v.dat");
    cout << "Nuclear Attraction Integral" << endl <<T << endl;
    cout << endl;
    // Forming Core Hamiltonian
    Eigen::Matrix<double, nao, nao> H = T + V;
    cout << "Core Hamiltonian" << endl << H << endl;
    cout << endl;

}

double hartree_fock::read_enuc(std::string int_path)
// Reading nuclear repulsion energy
// Expecting 'enuc.dat' at int_path
{
    std::string file_path = int_path + "/enuc.dat";
    std::ifstream enuc_file(file_path);
    //enuc = 0;
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
    int nlines = count_lines(filename);
    std::ifstream int_input(filename);
    if(int_input.is_open()) {
        int a, b; double temp;
        for (int i=0; i < nlines; i++)
        {
            int_input >> a >> b >> temp;
            one_e_matrix(a-1, b-1) = temp;
            one_e_matrix(b-1, a-1) = temp;
        }
    }
    else
    {
        cout << "Unable to open one-electron integral file" << endl;
    }

    return one_e_matrix;

}
//Eigen::MatrixXd hartree_fock::read_2e_ints(std::string int_path, std::string int_file)
//{
//    Eigen::MatrixXd two_e_matrix = Eigen::MatrixXd::Zero();
//}


// Should move this function to a general header file
int hartree_fock::count_lines(std::string filename)
{
    int nlines = 0;
    std::string line;
    std::ifstream input (filename);
    while(std::getline(input, line))
        ++ nlines;
    return nlines;
}