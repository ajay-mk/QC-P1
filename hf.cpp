// Hartree-Fock algorithm for practice
// Integrals are read from external files in this version
#include "hf.h"
#include "molecule.h"
#include "general.h"
#include <string>
#include <fstream>
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

// SCF Loop Parameters
const int max_iter = 10;
const double conv = 1e-10;


// Main Hartree-Fock Function
hartree_fock::hartree_fock(std::string input, const std::string int_path){
    Molecule mol(input);
    mol.print_info();
    const int nocc = mol.nelectrons()/2;
    cout << "Number of occupied states: " << nocc << endl;
    mol.print_geometry();

    // Reading Nuclear Repulsion
    double enuc  = read_enuc(int_path);
    cout << "----Nuclear repulsion energy----" << endl << enuc << " Eh" <<  endl << endl;

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
    auto Hcore = T + V;
    cout << "----Core Hamiltonian----" << endl << Hcore << endl;
    cout << endl;

    // Reading 2-e integral
    Eigen::MatrixXd two_e = read_2e_ints(int_path, "/eri.dat");

    // INPUT READING COMPLETE

    // STARTING CALCULATIONS

    // Diagonalizing S
    hartree_fock::diag_results S_Diag;
    S_Diag = Diag_M(S, S_Diag);
    auto S_evecs = S_Diag.evecs;
    auto S_evals = S_Diag.evals;

    // Building and Printing Symmetric Orthogonolization Matrix
    for (int i=0; i < S.rows(); i++)
    {
        S_evals(i, i) = pow(S_evals(i, i), -0.5); // Square root of S_evals
    }
    auto SOM = S_evecs * S_evals * S_evecs.transpose();
    cout << "----Symmetric Orthogonalization Matrix----" << endl;
    cout << SOM << endl;
    cout << endl;

    // Building and Printing Initial Fock Matrix - Using Core Hamiltonian
    auto F_initial = SOM.transpose() * Hcore * SOM ;
    cout << "----Initial Fock Matrix----" << endl;
    cout << F_initial << endl;

    //Diagonalizing Initial Fock Matrix
    hartree_fock::diag_results F_Diag;
    F_Diag = Diag(F_initial, F_Diag);
    auto C0_ = F_Diag.evecs;
    auto eps = F_Diag.evals;

    // Transformation of Eigenvectors into original AO basis
    auto C = SOM * C0_;
    cout << endl;
    cout << "----Initial Coefficient Matrix----" << endl;
    cout << C << endl;

    // Building and Printing Initial Density Matrix
    Eigen::MatrixXd D = build_density(C, nocc);
    cout << endl;
    cout << "----Initial Density Matrix----" << endl;
    cout << D << endl;

    // Computing Initial SCF energy
    double e_hf = scf_energy(D, Hcore, F_initial);

    double e_tot = e_hf + enuc ;
    cout << endl;
    cout << "Initial electronic Energy: " << e_hf << " Eh" << endl;
    cout << "Initial SCF Energy: " << e_tot << " Eh" << endl;

//    do {
//        auto D_old = D_n;
//        e_last = e_hf;
//        ++iter;
//        auto F_n = fock_build(Hcore, D, two_e);
//
//        F_n = SOM.transpose() * F_n * SOM;
//
//        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver2(F_n);
//        auto C_n = solver2.eigenvectors();
//        auto eps_n = solver2.eigenvalues();
//
//        auto C = SOM * C_n;
//
//        auto D = build_density(C, nocc);
//
//        e_hf = scf_energy(D, Hcore, F_n);
//        cout << "Iter " << iter << endl;
//        cout << e_hf << endl;
//        e_diff = e_hf - e_last;
//        cout << "E_diff: " << e_diff << endl;
//        D_old = D;
//        e_last = e_hf;
//
//    } while ((iter < max_iter) || (fabs(e_diff) > conv));

    // SCF Loop
    int iter = 0;
    double e_diff;
    do {
        auto ehf_last = e_hf;
        auto D_last = D;
        ++iter;
        // New Fock Matrix
        auto F = fock_build(Hcore, D_last, two_e);

        // Orthogonalize
        F = SOM.transpose() * F * SOM;

        // Diagonalize
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver2(F);
        auto C_n = solver2.eigenvectors();
        auto eps_n = solver2.eigenvalues();

        // Transform
        auto C = SOM * C_n;

        // New Density Matrix
        auto D = build_density(C, nocc);

        // New HF energy
        e_hf = scf_energy(D, Hcore, F);

        // Energy difference
        e_diff = e_hf - ehf_last;

        if (iter == 1)
            std::cout <<
                      "\n\n Iter        E(elec)              E(tot)               Delta(E)\n";
        printf(" %02d %20.12f %20.12f %20.12f\n", iter, e_hf, e_hf + enuc, e_diff);

    } while (((fabs(e_diff) > conv)) && (iter < max_iter));


}
// Function Definitions

double hartree_fock::read_enuc(std::string int_path)
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
hartree_fock::diag_results hartree_fock::Diag(Eigen::MatrixXd M, diag_results) {
    diag_results results;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(M);
    results.evecs = solver.eigenvectors();
    results.evals = solver.eigenvalues();
    return results;
}


hartree_fock::diag_results hartree_fock::Diag_M(Eigen::MatrixXd M, diag_results) {
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

Eigen::MatrixXd hartree_fock::build_density(Eigen::MatrixXd Coeff, int nocc) {
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

double hartree_fock::scf_energy(Eigen::MatrixXd D, Eigen::MatrixXd H,  Eigen::MatrixXd F)
{
    double e_hf = 0;
    for(int i=0; i < nao; i++){
        for(int j=0; j < nao; j++){
            e_hf += D(i, j) * (H(i, j) + F (i, j));
        }
    }
    return e_hf;
}

Eigen::MatrixXd hartree_fock::fock_build(Eigen::MatrixXd H, Eigen::MatrixXd D, Eigen::MatrixXd two_e)
{
    auto F = H;
    for (int i = 0; i < nao; i++) {
        for (int j = 0; j < nao; j++) {
            for (int k = 0; k < nao; k++) {
                for (int l = 0; l < nao; l++) {
                    int ij = CompoundIndex(i, j);
                    int kl = CompoundIndex(k, l);
                    int ik = CompoundIndex(i, k);
                    int jl = CompoundIndex(j, l);
                    F(i, j) += D(k, l) * (2 * two_e(CompoundIndex(ij, kl)) - two_e(CompoundIndex(ik, jl)));
                }
            }
        }
    }
    return F;
}