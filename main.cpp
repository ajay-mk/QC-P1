#include "hf.h"
#include "molecule.h"
#include <iostream>
#include <iomanip>

// Hartree-Fock algorithm for practice
// Integrals are read from external files in this version

using std::cout, std::endl;

int main(int argc, char *argv[]) {
  std::string filename = argv[1]; // Reading the filename from the argument
  std::string int_path = argv[2];
  cout << std::setprecision(12);
  cout << "Input file: " << filename << endl; // Printing the filename
  cout << "Integral files are in: " << int_path << endl;

  Molecule mol(filename);
  mol.print_info();
  const int nocc = mol.nelectrons() / 2;
  cout << "Number of occupied states: " << nocc << endl;
  mol.print_geometry();

  // Hartree-Fock
  // Since we are considering water molecule in STO-3G basis
  const int nao = 7;

  // Reading Nuclear Repulsion
  double enuc = read_enuc(int_path);
  cout << "\tNuclear repulsion energy" << endl << enuc << " Eh" << endl << endl;

  //Note: All the integrals are read in AO basis

  // Reading 1e integrals
  Eigen::Matrix<double, nao, nao> S = read_1e_ints(int_path, "/overlap.dat");
  cout << "\tOverlap Integrals" << endl << S << endl;
  cout << endl;
  Eigen::Matrix<double, nao, nao> T = read_1e_ints(int_path, "/ke.dat");
  cout << "\tKinetic Energy Integrals" << endl << T << endl;
  cout << endl;
  Eigen::Matrix<double, nao, nao> V = read_1e_ints(int_path, "/v.dat");
  cout << "\tNuclear Attraction Integrals" << endl << V << endl;
  cout << endl;
  // Forming Core Hamiltonian
  auto Hcore_ao = T + V;
  cout << "\tCore Hamiltonian" << endl << Hcore_ao << endl;
  cout << endl;

  // Reading 2-e integral
  Eigen::MatrixXd two_e = read_2e_ints(int_path, "/eri.dat");

  // Diagonalizing S
  diag_results S_Diag;
  S_Diag = Diag_M(S);
  auto S_evecs = S_Diag.evecs;
  auto S_evals = S_Diag.evals;

  // Building and Printing Symmetric Orthogonolization Matrix
  for (int i = 0; i < S.rows(); i++) {
    S_evals(i, i) = pow(S_evals(i, i), -0.5); // Square root of S_evals
  }

  auto SOM_oao = S_evecs * S_evals * S_evecs.transpose();
  cout << "\tSymmetric Orthogonalization Matrix" << endl;
  cout << SOM_oao << endl;
  cout << endl;
  // SOM is in Orthogonal AO (OAO) basis

  // Building and Printing Initial Fock Matrix - Using Core Hamiltonian
  auto F_oao = SOM_oao.transpose() * Hcore_ao * SOM_oao;
  cout << "\tInitial (Guess) Fock Matrix" << endl;
  cout << F_oao << endl;

  // Diagonalizing Initial Fock Matrix
  diag_results F_Diag;
  F_Diag = Diag(F_oao);
  auto C_oao = F_Diag.evecs;
  auto eps = F_Diag.evals;

  // Transformation of Eigenvectors into original AO basis
  auto C_mo = SOM_oao * C_oao;
  cout << endl;
  cout << "\tInitial Coefficient Matrix" << endl;
  cout << C_mo << endl;

  // Building and Printing Initial Density Matrix
  Eigen::MatrixXd D_ao = build_density(C_mo, nocc);
  cout << endl;
  cout << "\tInitial (Guess) Density Matrix" << endl;
  cout << D_ao << endl;

  // Computing Initial SCF energy
  double e_hf0 = scf_energy(D_ao, Hcore_ao, Hcore_ao);

  cout << endl;
  cout << "From Initial Guess: " << endl;
  cout << "Initial electronic Energy: " << e_hf0 << " Eh" << endl;
  cout << "Initial Total Energy: " << e_hf0 + enuc << " Eh" << endl << endl;

  // SCF Loop Parameters
  const int max_iter = 100;
  const double conv = 1e-12;

  // SCF Loop

      int iter = 0;
      double e_diff;
      double e_hf = e_hf0;
      do {

          ++iter;
          //cout << iter << endl;
          auto ehf_last = e_hf;

          // New Fock Matrix
          auto F_ao = fock_build(Hcore_ao, D_ao, two_e);
          // Orthogonalize
          auto F_oao = SOM_oao.transpose() * F_ao * SOM_oao;
          if(iter ==1)
              cout << "Fock Matrix" << endl << F_oao << endl << endl;

          // Diagonalize
          diag_results Fn_diag;
          Fn_diag = Diag(F_oao);
          C_oao = Fn_diag.evecs;
          auto eps_n = Fn_diag.evals;

          // Transform
          auto C_mo = SOM_oao * C_oao;

          // New Density Matrix
          D_ao = build_density(C_mo, nocc);
          //cout << "Testing Density" << endl << D << endl << endl;

          // New HF energy
          e_hf = scf_energy(D_ao, Hcore_ao, F_ao);

          // Energy difference
          e_diff = abs(e_hf - ehf_last);


          if (iter == 1)
              std::cout <<
                        "\n\n Iter        E(elec)              E(tot)                Delta(E)\n";
          printf(" %02d %20.12f %20.12f %20.12f\n", iter, e_hf, e_hf + enuc,e_diff);
      } while ((fabs(e_diff) > conv) && (iter < max_iter));

      cout << endl
           << "Hartree Fock Energy = " << e_hf + enuc << " Eh" << endl;

  return 0;
}