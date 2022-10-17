#include "molecule.h"
#include <fstream>
#include <cmath>
#include "mass.h"

using std::cout, std::endl;

// Defining molecule class
Molecule::Molecule(const std::string filename){
    //Opening input file/geometry
    std::ifstream input (filename);
    if(input.is_open()){
        input >> natoms;
    }
    else
    {
        cout << "Unable to open " << filename << endl;
    }

    // Reading Geometry
    atoms.resize(natoms);
    for(int i=0; i < natoms; i++) {
        // Storing coordinates into vector
        input >> atoms[i].Z >> atoms[i].x >> atoms[i].y >> atoms[i].z;
        // Calculating number of electrons
        nelectrons += atoms[i].Z;
        // Calculating molecular mass
        mol_mass += mass[atoms[i].Z];
    }
}

void Molecule::print_geometry()
{
    cout << endl;
    for(int i=0; i < natoms; i++)
    {
        printf("%d\t%4.10f\t%4.10f\t%4.10f\n", atoms[i].Z, atoms[i].x, atoms[i].y, atoms[i].z);
    }
}

void Molecule::print_info()
{
    cout << endl;
    cout << "Number of atoms: " << natoms<< endl;
    cout << "Number of electrons: " << nelectrons << endl;
    cout << "Molecular Mass: " << mol_mass << " amu" << endl;
}

double Molecule::bond_length(int i, int j) {
    double bl = sqrt(pow((atoms[i].x-atoms[j].x),2)
                     + pow((atoms[i].y-atoms[j].y),2)
                     + pow((atoms[i].z-atoms[j].z),2));
    return bl;
}

void Molecule::bond_matrix()
{
    cout << "Bond Lengths" << endl;
    printf("i\tj\tBond Length");
    bonds.resize(natoms, natoms);
    for(int i=0; i < natoms; i++){
        for(int j = 0; j < natoms; j++){
            bonds(i, j) = Molecule::bond_length(i, j);
            if(bonds(i, j) < 4){
            // There are two ways two print this, print atomic numbers or indices
            //printf("%d\t%d\t%4.5f\n", atoms[i].Z, atoms[j].Z, bonds(i,j));
            printf("%d\t%d\t%4.5f\n", i, j, bonds(i,j));
            }
        }
    }
    //cout << "Bond Length Matrix" << endl;
    //cout << bonds << endl;
}

void Molecule::translate(double x, double y, double z)
{
    printf("Translating molecule by (%4.3f, %4.3f, %4.3f)\n", x, y, z);
    for(int i = 0; i < natoms; i++)
    {
        atoms[i].x += x;
        atoms[i].y += y;
        atoms[i].z += z;
    }
    cout << "Printing updated geometry:" << endl;
    Molecule::print_geometry();
}
double Molecule::bond_angle(int i, int j, int k)
// cos phi(ijk) = u(ji).u(jk), where u is an unit vector
// u = (-(xi - xj)/Rij, -(yi-yj/Rij, -(zi-zj)/Rij)) --> See the unit_vector function
{
    double dot_pdt = unit_vector(j,i).dot(unit_vector(j,k));
    double angle_radian = acos(dot_pdt);
    double angle_degree = angle_radian * 180/M_PI;
    return angle_degree;
}

Eigen::Vector3d Molecule::unit_vector(int i, int j)
{
    Eigen::Vector3d uv = { (-(atoms[i].x-atoms[j].x)/bond_length(i, j)),
                           (-(atoms[i].y-atoms[j].y)/ bond_length(i, j)),
                           (-(atoms[i].z-atoms[j].z)/ bond_length(i, j))};
    return uv;
}

void Molecule::ba_matrix()
{
    cout << "Printing Bond Angles" << endl;
    printf("i\tj\tk\tphi\n");
    for(int i=0; i < natoms; i++){
        for(int j=0; j < i ; j++){
            for(int k=0; k < j; k++)
            {
                // There are two ways two print this, print atomic numbers or indices
                printf("%d\t%d\t%d\t%4.5f\n", i, j, k, bond_angle(i, j, k));
                //printf("%d\t%d\t%d\t%4.5f\n", atoms[i].Z, atoms[j].Z, atoms[k].Z, bond_angle(i, j, k));
            }
        }
    }

}

Eigen::Vector3d Molecule::find_com()
// Computing centre of mass coordinates
// Xcm = Sum(mixi)/total mass
{
    double Nr_x, Nr_y, Nr_z;
    for(int i =0; i < natoms; i++)
    {
        Nr_x =+ atoms[i].x * mass[atoms[i].Z];
        Nr_y =+ atoms[i].y * mass[atoms[i].Z];
        Nr_z =+ atoms[i].z * mass[atoms[i].Z];
    }
    com = {Nr_x/mol_mass, Nr_y/mol_mass, Nr_z/mol_mass};
    printf("Center of mass coordinates: (%4.5f, %4.5f, %4.5f)\n", com(0), com(1), com(2));
    return com;
}
