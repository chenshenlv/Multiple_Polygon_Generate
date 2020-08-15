#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <cstdlib>  // For srand() and rand()

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <boost/foreach.hpp>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <boost/function_output_iterator.hpp>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/intersections.h>
#include <CGAL/Polygon_mesh_processing/smooth_mesh.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>



typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                    Point;
typedef Kernel::Vector_3                   Vector;
typedef CGAL::Surface_mesh<Point>          Surface_mesh;
typedef CGAL::Aff_transformation_3<Kernel> Aff_transformation_3;
typedef CGAL::Point_2<Kernel> Point_2;
typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef CGAL::Aff_transformation_3<Kernel> Aff_transformation_3;
typedef Kernel::RT                                   RT_num;
typedef Kernel::Segment_2 Segment_2;
typedef Kernel::Intersect_2 Intersect_2;


typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Surface_mesh>::edge_descriptor     edge_descriptor;

typedef boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Surface_mesh>::face_descriptor   face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;
using namespace std;
extern Surface_mesh mesh;
struct halfedge2edge
{
  halfedge2edge(const Surface_mesh& m, std::vector<edge_descriptor>& edges)
    : m_mesh(m), m_edges(edges)
  {}
  void operator()(const halfedge_descriptor& h) const
  {
    m_edges.push_back(edge(h, m_mesh));
  }
  const Surface_mesh& m_mesh;
  std::vector<edge_descriptor>& m_edges;
};

void ISO_MESH(int mesh_no, double target_edge_length, Surface_mesh &mesh, int EDGE_POLY);
void transform(int mesh_no,double radius,Surface_mesh &mesh);

#endif