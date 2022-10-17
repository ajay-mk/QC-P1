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
        //mol.bond_matrix();
        //mol.translate(2,2,2);
        //mol.ba_matrix();
        mol.find_com();
        cout << endl;



    return 0;
}
}