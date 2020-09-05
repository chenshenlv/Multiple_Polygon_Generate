#include <stdio.h>
#include <stdlib.h>  
#include <iostream>
#include "polygon.h"
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <getopt.h>
#include <new>

const double C = 343.21;

int main(int argc, char** argv)
{ 
  int c;
  int option_index = 0;
  // string filename;
  string outfiles = "merged_mesh.off";
  // filename = "trans_mesh1.off";
  // stringstream a;
  // a << mesh_no;
  
  std::ofstream out;
  std::ifstream input;
  int total_scene_no = 5; // no. of result obj files
  int max_mesh_num_in_scene = 3;
  int edge_lower = 3; // no. of polygon edge low range
  int edge_upper = 7; // no. of polygon edge upper range
  double RADIUS = 0.5; // radius for polygon excircle
  double HIGHT = 1.0; // hight of extruded polygon
  double freq_max = 2000; // max frequency of sound
  double target_edge_length = 0.02; // define target_edge_length of mesh 
  string output_file_path = "mesh/";
  int EDGE_POLY = 9;
  

  // static struct option long_options[] = {
  //   {"total_scene_no",  required_argument, 0,  0 },
  //   {"max_mesh_num_in_scene",   required_argument, 0,  1 },
  //   {"edge_lower ",   required_argument, 0,  2 },
  //   {"edge_upper",    required_argument, 0,  3 },
  //   {"radius",        required_argument, 0,  4 },
  //   {"height",        required_argument, 0,  5 },
  //   {"freq_max",      required_argument, 0,  6 },
  //   {"output_path",   required_argument, 0,  7 },
  //   {0,             0,                 0,  0 }
  // };
  
  // while(1){
  //   c = getopt_long(argc, argv, "", long_options, &option_index);
  //   // std::cout<< "c"<<c<<std::endl;
  //   if(c == -1)
  //     break;

  //   switch(c){

  //     case 0:
  //           total_scene_no = atoi(optarg);
  //           // std::cout<<"total_scene_no = "<< total_scene_no << std::endl;
  //           break;
  //     case 1:
  //           max_mesh_num_in_scene = atoi(optarg);
  //           // std::cout<<"max_mesh_num_in_scene = "<< max_mesh_num_in_scene << std::endl;
  //           break;
  //     case 2:
  //           edge_lower = atoi(optarg);
  //           // std::cout<<"edge_lower = "<< edge_lower << std::endl;
  //           break;
  //     case 3:
  //           edge_upper = atoi(optarg);
  //           // std::cout<<"edge_upper = "<< edge_upper << std::endl;
  //           break;
  //     case 4:
  //           RADIUS = atof(optarg);
  //           // std::cout<<"radius = "<< RADIUS << std::endl;
  //           break;
  //     case 5:
  //           HIGHT = atof(optarg);
  //           // std::cout<<"extrude_height = "<< HIGHT << std::endl;
  //           break;
  //     case 6:
  //           freq_max = atof(optarg);
  //           target_edge_length = C/(7*freq_max);
  //           // std::cout<<"freq_max = "<< freq_max << std::endl;
  //           // std::cout<<"target_edge_length = "<< target_edge_length << std::endl;
  //           break;
  //     case 7:
  //           output_file_path = optarg;
  //           // std::cout<<"output_file_path = "<< output_file_path<< std::endl;
  //           break;
  //     default:
  //           std::cerr<<"illegal arguments inputed"<<std::endl;
  //           // std::cout<< "c"<<c<<std::endl;
  //   }
  // }

  // if (argc < sizeof(long_options)/sizeof(long_options[0])){
  //   std::cerr<<"argments input less than required, left arguments set as default"<<std::endl;
  //   std::cout<<"\ntotal_scene_no = "<< total_scene_no <<"\nedge_lower = "<< edge_lower <<"\nedge_upper = "<< edge_upper 
  //   <<"\nradius = "<< RADIUS<<"\nextrude_height = "<< HIGHT<<"\nfreq_max = "<< freq_max
  //   <<"\ntarget_edge_length = "<< target_edge_length<<std::endl;

  // }
  // if (edge_lower>edge_upper){
  //   edge_lower = edge_upper;
  //   std::cerr<<"edge low limit is bigger than up limit, set to up limit"<<std::endl;
  // }
  // // if (edge_no<=0){
  // //   std::cerr<<"edge range should be positive."<<std::endl;
  // //   return 1;
  // // }
  // if ( total_scene_no<=0 || HIGHT <= 0 || RADIUS <= 0 || freq_max <= 0 || edge_lower<=0 || edge_upper<=0){
  //   std::cerr<<"wrong parameters entered, zero or negtive being inputed." <<std::endl;
  //   return 1; 
  // }
  srand(time(0));
  int poly_no = 0;
  int side_point = 4;
  int side_up_bott = 2;
  for (int i = 0; i < total_scene_no; i++) //iterate files 
  {
    GENPOLY:
    Surface_mesh merged_mesh;
    for (int mesh_generate_iter = 0; mesh_generate_iter < max_mesh_num_in_scene; mesh_generate_iter++) //iterate generating meshes in each scene
    {
      
      int EDGE_POLY = (rand() % (edge_upper - edge_lower + 1)) + edge_lower; // No. of edges of polygon
      // std::cout << "edge generate "<<EDGE_POLY<<std::endl;
      Pt* vertex = new (nothrow) Pt[2*EDGE_POLY]();
      if (vertex == nullptr)
      std::cout << "Error: memory could not be allocated"<<std::endl;

      int** face = new int*[EDGE_POLY]();
      for(int j = 0; j < EDGE_POLY; ++j) {
          face[j] = new (nothrow) int[side_point]();
      }
      if (face == nullptr)
      std::cout << "Error: memory could not be allocated"<<std::endl;

      int** top_bott = new int*[EDGE_POLY]();
      for(int j = 0; j < EDGE_POLY; ++j) {
          top_bott[j] = new (nothrow) int[side_up_bott]();
      }
      if (top_bott == nullptr)
      std::cout << "Error: memory could not be allocated"<<std::endl;
      
      // int EDGE_POLY = 9;
      // std::cout << "generate polygon"<<std::endl;
      
      
      polygon( EDGE_POLY, RADIUS, HIGHT, max_mesh_num_in_scene, vertex, face, top_bott, poly_no);
      Surface_mesh mesh;
      output( 1, EDGE_POLY, vertex, face, top_bott, mesh );
      if (CGAL::Polygon_mesh_processing::is_outward_oriented(mesh))
      {
        CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh);
        std::cout<< "mesh normal is outward: " << CGAL::Polygon_mesh_processing::is_outward_oriented(mesh)<<std::endl;
        std::ofstream log;
        log.open ("log.txt", ios::out | ios::app ); 
        log<<" file no: "<<i<<std::endl;
        log<< "mesh normal is outward: "<< CGAL::Polygon_mesh_processing::is_outward_oriented(mesh)<<std::endl;
        log.close();
      }
      std::cout<< "mesh normal is outward: " << CGAL::Polygon_mesh_processing::is_outward_oriented(mesh)<<std::endl;
      
      // std::cout << "polygon done"<<std::endl;
      // ISO_MESH(1, target_edge_length);
      // std::cout<<"mesh gene:"<<mesh_generate_iter<<std::endl;
      for(int transform_iter = 0; transform_iter<50; transform_iter++)//iterate transforming single meshes
      {
        
        
        transform(1,RADIUS,mesh);
        // input.open(filename);
        // if (!input || !(input >> mesh))
        // {
        //   std::cerr << "First mesh is not a valid off file." << std::endl;
        
        //   return 1;
        // }
        // input.close();

        if (!CGAL::Polygon_mesh_processing::do_intersect(merged_mesh,mesh)) 
        { 
          ISO_MESH(i, target_edge_length,mesh,EDGE_POLY);
          CGAL::copy_face_graph(mesh, merged_mesh);
          // std::remove(filename.c_str());
          // std::remove("iso_mesh1.off");
          // std::cout<<transform_iter<<std::endl;
          break;
        }
        if (CGAL::Polygon_mesh_processing::do_intersect(merged_mesh,mesh) && transform_iter == 20)
        {
          delete [] vertex;

          for(int j = 0; j < EDGE_POLY; j++) {
            delete [] face[j];
          }
          delete [] face;

          for(int j = 0; j < EDGE_POLY; j++) {
            delete [] top_bott[j];
          }
          delete [] top_bott;
          std::cout<<"polygon too big, go back to GENPOLY"<<std::endl;
          goto GENPOLY;
        }
      }

      delete [] vertex;

      for(int j = 0; j < EDGE_POLY; j++) {
        delete [] face[j];
      }
      delete [] face;

      for(int j = 0; j < EDGE_POLY; j++) {
        delete [] top_bott[j];
      }
      delete [] top_bott;
    }

    Surface_mesh ground;
    double length = 1, width = 1, height = 0.1;
    create_ground(length, width, height, target_edge_length, ground);

    // if (!CGAL::Polygon_mesh_processing::do_intersect(merged_mesh,ground)) 
    // { 
    //   CGAL::copy_face_graph(ground, merged_mesh);
          
    // }
    Surface_mesh result;
    if (CGAL::Polygon_mesh_processing::do_intersect(merged_mesh,ground)) 
    { 
      
      bool valid_union = PMP::corefine_and_compute_intersection(merged_mesh,ground,result);
      
          
    }
    out.open("merged_mesh.off"); 
    out << result;
    out.close();
    off2obj(i,output_file_path);
    poly_no+=1;
    if (poly_no % 1000 == 0){std::cout<< "generated"<<poly_no<<"/"<<total_scene_no<< std::endl;}
      
    // std::remove(outfiles.c_str());
  }
  
  

}