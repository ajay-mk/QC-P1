#include <iostream>
#include "molecule.h"

using std::cout, std::endl;

int main(int argc , char* argv[]) {
    {
        std::string filename = argv[1]; // Reading the filename from the argument

        cout << "Input file: " << filename << endl; // Printing the filename

        Molecule mol(filename);
        mol.print_info();
        cout << endl;
        cout << "Printing Molecular Geometry:" << endl;
        mol.print_geometry();
        cout << endl;

        // Hartree-Fock

        //Inertia Tensor calculation needs debugging
        //Something is wrong while translating molecule to CoM --> Check

    return 0;
}
}
