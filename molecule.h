#include <iostream>
#include <cstdio>
#include <string>
#include <cmath>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>


// A Struct for information regarding geometry
// Contains atomic numbers and coordinates
struct Geom {
    int Z;
    double x, y, z;
};

// Class for managing molecular geometry and information
class Molecule {
public:
    Molecule(std::string);

    void print_info();
    void print_geometry();
    double bond_length(int a, int b);
    void distance_matrix();
    void translate (double x, double y, double z);
    //void rotate (double phi);
    double bond_angle (int a, int b, int c);
    void ba_matrix();
//    void oop_angles(int a, int b, int c);
//    void torsion (int a, int b, int c, int d);
    Eigen::Vector3d unit_vector(int a, int b);
    Eigen::Vector3d find_com();
    Eigen::MatrixXd compute_inertia_tensor();
    void moment_of_inertia();

private:
    int natoms;
    int total_charge;
    double mol_mass;
    double nelectrons;
    std::vector<Geom> atoms;
    Eigen::MatrixXd bonds;
    Eigen::Vector3d com;
    Eigen::MatrixXd inertia_tensor;
};

