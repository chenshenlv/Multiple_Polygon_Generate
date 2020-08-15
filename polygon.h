#ifndef POLYGON_H
#define POLYGON_H

#include"mesh.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
//#include <fstream>
#include <unistd.h>
#include <ctime>    // For time()
#include <cstdlib>  // For srand() and rand()
#include <cstddef>
//#include <cstring>
#include <stdio.h>     /* for printf */
#include <stdlib.h>    /* for exit */
#include <getopt.h>
#include <new>
#include <vector>

#include <cstring>
#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/IO/File_writer_wavefront.h>
#include <CGAL/IO/generic_copy_OFF.h>

const double PI = 3.141592653589793;
// using namespace std;

struct Pt{ double x, y, z; };

Pt operator +( Pt p, Pt q );
Pt operator -( Pt p, Pt q );
Pt operator *( double r, Pt p );
Pt dot_product(Pt p, Pt q);
double dot( Pt p, Pt q );
double abs( Pt p );
Pt cross( Pt p, Pt q );
Pt rotX( Pt p, double angle );       // Rotate about x axis
Pt rotZ( Pt p, double angle );       // Rotate about z axis
double distance_pol_cir(Pt p,  double radius); //distance between vertex to outter circile
Pt reflectInPlane( Pt p, Pt norm, Pt ref );    // Reflect in plane defined by normal vector and reference point
void polygon( int EDGE_POLY, double RADIUS, double HIGHT, int mesh_num, Pt *vertex, int **face, int **top_bott, int poly_no  );
void output( int mesh_no, int EDGE_POLY, Pt* vertex, int** face, int** top_bott, Surface_mesh &mesh );
void off2obj(int mesh_no,string output_file_path);   //convert off file to obj file
void create_ground ( double length, double width, double height, double target_edge_length, Surface_mesh &mesh);

#endif