#include <iostream>
#include "molecule.h"
#include "hf.h"
#include "general.h"

using std::cout, std::endl;

int main(int argc , char* argv[]) {
    {
        std::string filename = argv[1]; // Reading the filename from the argument
        std::string int_path = argv[2];
        cout << "Input file: " << filename << endl; // Printing the filename
        cout << "Integral files are in: " << int_path << endl;

        // Hartree-Fock
        hartree_fock calc_hf (filename, int_path);

    return 0;
}
}
