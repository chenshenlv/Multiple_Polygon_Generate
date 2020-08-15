#include "mesh.h"

void ISO_MESH(int mesh_no, double target_edge_length,Surface_mesh &mesh,int EDGE_POLY)
{   
    double edge_length_max, edge_length, area;
    // string filename;
    // string outfiles;
    // stringstream a;
    // a << mesh_no;
    // filename = "mesh" + a.str() + ".off";
    // outfiles = "iso_mesh" + a.str() + ".off";
    edge_length_max = 0;
    // const char* outfilename = (argc > 2) ? argv[2] : "P_tri.off";
    // std::ifstream input(filename);
    // Surface_mesh mesh;
    

    if ( mesh.is_empty())//
    {
    std::cerr << "Not a valid off file." << std::endl;

    }


    // CGAL::Polygon_mesh_processing::triangulate_faces(mesh);

    // Confirm that all faces are triangles.
    BOOST_FOREACH(boost::graph_traits<Surface_mesh>::face_descriptor fit, faces(mesh))
      if (next(next(halfedge(fit, mesh), mesh), mesh)
          !=   prev(halfedge(fit, mesh), mesh))
        std::cerr << "Error: non-triangular face left in mesh." << std::endl;

    //Compute edge_length
    BOOST_FOREACH(boost::graph_traits<Surface_mesh>::face_descriptor fit, faces(mesh))
      edge_length=CGAL::Polygon_mesh_processing::edge_length(halfedge(fit, mesh),mesh);
      if (edge_length>edge_length_max) edge_length_max = edge_length;
     
    unsigned int nb_iter = 5;

    //std::cout << "Split border...";
    
    CGAL::Timer t;
    t.start();
    std::vector<edge_descriptor> border;
    PMP::border_halfedges(faces(mesh),
    mesh,
    boost::make_function_output_iterator(halfedge2edge(mesh, border)));
    PMP::split_long_edges(border, target_edge_length, mesh);

    //std::cout << "done." << std::endl;

    // std::cout << "Start remeshing of " << filename
    // << " (" << num_faces(mesh) << " faces)..." << std::endl;
    // Constrain edges with a dihedral angle over 60°
    typedef boost::property_map<Surface_mesh, CGAL::edge_is_feature_t>::type EIFMap;
    EIFMap eif = get(CGAL::edge_is_feature, mesh);
    PMP::detect_sharp_edges(mesh, 10, eif);
    int sharp_counter = 0;
    for(edge_descriptor e : edges(mesh)){
      if(get(eif, e))
        ++sharp_counter;
    }
    int chk_ite=0;
    // std::cout << sharp_counter << " sharp edges" << std::endl;
    while (edge_length_max>target_edge_length*1.1){
      edge_length_max=0;

      PMP::isotropic_remeshing(
      faces(mesh),
      target_edge_length,
      mesh,
      PMP::parameters::number_of_iterations(nb_iter)
      .edge_is_constrained_map(eif)
      //.protect_constraints(true)//i.e. protect border, here
      );

      //Compute edge_length
      BOOST_FOREACH(boost::graph_traits<Surface_mesh>::face_descriptor fit, faces(mesh))
      edge_length=CGAL::Polygon_mesh_processing::edge_length(halfedge(fit, mesh),mesh);
      if (edge_length>edge_length_max) edge_length_max = edge_length;
      chk_ite++;
      if (chk_ite>5) break;

      }

    // // Constrain edges with a dihedral angle over 60°
    // typedef boost::property_map<Surface_mesh, CGAL::edge_is_feature_t>::type EIFMap;
    // EIFMap eif = get(CGAL::edge_is_feature, mesh);
    // PMP::detect_sharp_edges(mesh, 60, eif);
    // int sharp_counter = 0;
    // for(edge_descriptor e : edges(mesh))
    //   if(get(eif, e))
    //     ++sharp_counter;
    // std::cout << sharp_counter << " sharp edges" << std::endl;
    // // Smooth with both angle and area criteria + Delaunay flips
    // PMP::smooth_mesh(mesh, PMP::parameters::number_of_iterations(10)
    //                                      .use_safety_constraints(false) // authorize all moves
    //                                      .edge_is_constrained_map(eif));

    
    // Compute surface area.
    area=PMP::area(faces(mesh),mesh,
    PMP::parameters::vertex_point_map(get(CGAL::vertex_point, mesh)));
    // std::cerr << "Done (" << t.time() << " s)" << " (" << mesh.number_of_faces() << " faces)"
    // <<" ( surface area:" << area << " )" << " ( max edge length " << edge_length_max << " )..."<< std::endl;
    // std::ofstream out(outfiles);
    // out << mesh;
    // std::cout << "Remeshing done.\n\n" << std::endl;

    // std::remove(filename.c_str());
   
}


void transform(int mesh_no,double radius,Surface_mesh &mesh)
{
	// string filename;
 //    string outfiles;
    // stringstream a;
    // a << mesh_no;
    // filename = "mesh" + a.str() + ".off";
    // outfiles = "trans_mesh" + a.str() + ".off";
    // std::ifstream input(filename);
    // Surface_mesh mesh;

    // if (!input || !(input >> mesh) || mesh.is_empty())//
    // {
    //   std::cerr << "Not a valid off file/trans." << std::endl;
      
    // }
    // input.close();
    Kernel::Iso_cuboid_3 c3 =CGAL::Polygon_mesh_processing::bbox (mesh); 
    double x_low = -c3[0][0]-radius;
    double x_up = radius - c3[1][0];
    double y_low = -c3[0][1]-radius;
    double y_up = radius - c3[2][1];
    double scale_x = ((double)rand() / double(RAND_MAX))*(x_up-x_low)+x_low; // scale polygon
    double scale_y = ((double)rand() / double(RAND_MAX))*(y_up-y_low)+y_low; // scale polygon
	Aff_transformation_3 translate(CGAL::TRANSLATION, Vector(scale_x, scale_y, 0));
    CGAL::Polygon_mesh_processing::transform (translate, mesh);  
    // std::cout<<"c3:"<<c3<<std::endl;
    // std::cout<<"c3[0]:"<<c3[0]<<std::endl;
    // std::cout<<"c3[1]:"<<c3[1]<<std::endl;
    // std::cout<<"c3[2]:"<<c3[2]<<std::endl;
    // std::cout<<"c3[3]:"<<c3[3]<<std::endl;
    // std::cout<<"c3[0][0]:"<<c3[0][0]<<std::endl;
    // std::cout<<"c3[0][1]:"<<c3[0][1]<<std::endl;
    // std::cout<<"c3[1][0]:"<<c3[1][0]<<std::endl;
    // std::cout<<"c3[1][1]:"<<c3[1][1]<<std::endl;

    // std::cout<<"x_rang["<<x_low<<","<<x_up<<"]"<<std::endl;
    // std::cout<<"scale_x"<<scale_x<<std::endl; 
    // std::cout<<"y_rang["<<y_low<<","<<y_up<<"]"<<std::endl;
    // std::cout<<"scale_y"<<scale_y<<std::endl; 
    // std::ofstream out(outfiles);
    // out << mesh;
    // std::remove(filename.c_str());
    
}