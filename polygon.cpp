#include "polygon.h"

Pt operator +( Pt p, Pt q ) { return { p.x + q.x, p.y + q.y, p.z + q.z }; }
Pt operator -( Pt p, Pt q ) { return { p.x - q.x, p.y - q.y, p.z - q.z }; }
Pt operator *( double r, Pt p ) { return { r * p.x, r * p.y, r * p.z }; }
Pt dot_product(Pt p, Pt q) { return {p.x * q.x, p.y * q.y, p.z * q.z};}
double dot( Pt p, Pt q ) { return p.x * q.x + p.y * q.y + p.z * q.z; }
double abs( Pt p ) { return sqrt( dot( p, p ) ); }
Pt cross( Pt p, Pt q ) { return { p.y * q.z - p.z * q.y, p.z * q.x - p.x * q.z, p.x * q.y - p.y * q.x }; }
Pt rotX( Pt p, double angle )       // Rotate about x axis
{
   double c = cos( angle ), s = sin( angle );
   return { p.x, p.y * c - p.z * s, p.y * s + p.z * c };
}


//================================================


Pt rotZ( Pt p, double angle )       // Rotate about z axis
{
   double c = cos( angle ), s = sin( angle );
   return { p.x * c - p.y * s, p.x * s + p.y * c, p.z };
}
//================================================
double distance_pol_cir(Pt p,  double radius) //distance between vertex to outter circile
{
  return (radius-abs(p));
}
//================================================

Pt reflectInPlane( Pt p, Pt norm, Pt ref )    // Reflect in plane defined by normal vector and reference point
{
   norm = ( 1.0 / abs( norm ) ) * norm;       // Make unit normal
   return { p - 2.0 * dot( p - ref, norm ) * norm };
}


void polygon( int EDGE_POLY, double RADIUS, double HIGHT, int mesh_num, Pt *vertex, int **face, int **top_bott, int poly_no  )
{
   LABEL1:

   vertex[0] = {RADIUS, 0, 0};
   Pt extrude = {0, 0, HIGHT};
   double scale_x, scale_y;
   double delta_x, delta_y;
   Pt scale;
   // int i = 1;
   double angle_left = 360;
   double angle_seed;
   double angle;
   // double x_range, y_range, x_range_min, y_range_min;
   double x_range_left, x_range_right, y_range_top, y_range_bot;
   
   double* angle_arr = new (nothrow) double[EDGE_POLY]();
      if (angle_arr == nullptr)
      std::cout << "Error: angle_arr memory could not be allocated";
   angle_arr[0]=0;
   // std::cout<< "Generate vertex..."<<std::endl; 
   // std::cout<<"vertex0"<<":"<<"("<<vertex[0].x<<" ,"<<vertex[0].y<<" ,"<<vertex[0].z<<")"<<std::endl;  
   
   for (int i = 1; i< EDGE_POLY; i++)
    {
       GEN_ANG:
               angle_seed = ((double)rand() / (double)RAND_MAX)*(0.99-0.027)+0.027;
               angle_arr[i] = 360*angle_seed;
               for (int j = 0; j < i ; j++)
               {
                if(abs(angle_arr[j]-angle_arr[i]) < 360*0.027) goto GEN_ANG;
               }

         
      // angle_seed = ((double)rand() / (double)RAND_MAX)*(0.99-0.2)+0.2;     
      // double angle = angle_seed*(angle_left/(EDGE_POLY-i)); //angle for generate vertics
      // //double angle = (rand() % ((int)(360*1.1/EDGE_POLY)-((int)(360/EDGE_POLY/3))+1))+((int)(360/EDGE_POLY/3));
      // angle_left = angle_left - angle;
      // vertex[i] = rotZ(vertex[i-1],angle * PI/180);
      // // std::cout<<"vertex"<<i<<":"<<"("<<vertex[i].x<<" ,"<<vertex[i].y<<" ,"<<vertex[i].z<<")"<<std::endl;
    }
    double temp;
    for (int j = 0; j< EDGE_POLY; j++)
    {
      for (int k = j+1; k < EDGE_POLY; k++)
      {
        if ( angle_arr[k]<angle_arr[j])
        {
          temp = angle_arr[j];
          angle_arr[j]=angle_arr[k];
          angle_arr[k]=temp;

        }

      }
    }

    for (int i = 1; i < EDGE_POLY; i++)
    {
      vertex[i] = rotZ(vertex[0],angle_arr[i] * PI/180);
      // std::cout<<"angle"<<i<<"="<<angle_arr[i]<<std::endl;
    }

    delete [] angle_arr;
     std::cout<<"gen_ang done"<<std::endl;     

   // std::cout<< "Generate Done"<<std::endl;
   // std::cout<< "Scale vertex..."<<std::endl; 
   double roate_angle = ((double)rand() / double(RAND_MAX))*(2*PI-0)+0; // Roatation angle
   for (int i = 0; i < EDGE_POLY; i++)
    {
      vertex[i] = rotZ(vertex[i],roate_angle);
      // std::cout<<"angle"<<i<<"="<<angle_arr[i]<<std::endl;
    }
   std::cout<<"roate_angle done"<<std::endl;
   
   // x_range_min = RADIUS;
   // y_range_min = RADIUS;
   // std::cout<<"scale_x"<<scale_x<<" ,"<<"scale_y"<<scale_y<<std::endl;
   for ( int i = 0; i < EDGE_POLY; i++)
   {
		// if ( scale_x > scale_min )  scale_min = scale_x;
		//scale_y = scale_x;
		// std::cout<<"before scale vertex"<<i<<":"<<"("<<vertex[i].x<<" ,"<<vertex[i].y<<" ,"<<vertex[i].z<<")"<<std::endl;
	 	scale_x = ((double)rand() / double(RAND_MAX))*(0.99-0.25)+0.25; // scale polygon
	 	scale_y = ((double)rand() / double(RAND_MAX))*(0.99-0.25)+0.25; // scale polygon
		scale={scale_x, scale_x, 0};
		vertex[i] = dot_product(vertex[i],scale);
		// std::cout<<"vertex"<<i<<":"<<"("<<vertex[i].x<<" ,"<<vertex[i].y<<" ,"<<vertex[i].z<<")"<<std::endl;
		// x_range = RADIUS - abs(vertex[i].x);
		// y_range = RADIUS - abs(vertex[i].y);
		// x_range_min = (x_range<x_range_min)?x_range:x_range_min;
		// y_range_min = (y_range<y_range_min)?y_range:y_range_min;
      
   }
   std::cout<<"scale done"<<std::endl;
   Polygon_2 bot_poly;
   bot_poly.clear();
   Pt* disper_point = new (nothrow) Pt[2500]();
   Segment_2* poly_seg = new (nothrow) Segment_2[EDGE_POLY]();
   if (poly_seg == nullptr)
   std::cout << "Error: poly_seg memory could not be allocated";
   
   // get 2D bot_poly and poly edge segment
   for (int j = 0; j < EDGE_POLY; j++)
   {

    bot_poly.push_back(Point_2(vertex[j].x,vertex[j].y));
    // double dot_len = sqrt(vertex[j].x*vertex[j].x+vertex[j].y*vertex[j].y);
    // std::cout << dot_len << endl;

    if (j+1<EDGE_POLY)
    {
      poly_seg[j] = Segment_2(Point_2(vertex[j].x,vertex[j].y),Point_2(vertex[j+1].x,vertex[j+1].y));
    }
    if(j+1==EDGE_POLY)
    {
      poly_seg[j] = Segment_2(Point_2(vertex[j].x,vertex[j].y),Point_2(vertex[0].x,vertex[0].y));
    }
    if (j+2<EDGE_POLY)
    {
        if (CGAL::collinear(Point_2(vertex[j].x,vertex[j].y),Point_2(vertex[j+1].x,vertex[j+1].y),Point_2(vertex[j+2].x,vertex[j+2].y)))
      {
        std::cout<<"collinear go back to LABEL1"<<std::endl;
        delete [] disper_point;
        delete [] poly_seg;
        goto LABEL1;
        
      }
      
    }
    if (j+2==EDGE_POLY)
    {
        if (CGAL::collinear(Point_2(vertex[j].x,vertex[j].y),Point_2(vertex[j-1].x,vertex[j-1].y),Point_2(vertex[j-2].x,vertex[j-2].y)))
      {
        std::cout<<"collinear go back to LABEL1"<<std::endl;
        delete [] disper_point;
        delete [] poly_seg;
        goto LABEL1;
        
      }
      
    }
    

   }
   std::cout<<"get 2D bot_poly and poly edge segment done"<<endl;

   //get bounding box
   // std::cout<<bot_poly.bbox()<<endl;
   // std::cout<<x_range_left<<endl;
   // std::cout<<x_range_right<<endl;
   // std::cout<<y_range_top<<endl;
   // std::cout<<y_range_bot<<endl;



   // std::cout<<"bot_poly done"<<std::endl;
   // CGAL::set_pretty_mode(std::cout);
   // std::cout << "created the polygon p:" << std::endl;
   // std::cout << bot_poly << std::endl;
   // std::cout << std::endl;
   // std::cout << Point_2(vertex[0].x,vertex[0].y) << std::endl;
   // std::cout << Point_2(vertex[1].x,vertex[1].y) << std::endl;
   // std::cout << Point_2(vertex[2].x,vertex[2].y) << std::endl;
   // std::cout << std::endl;
   // determine some properties of the polygon
   bool IsSimple    = bot_poly.is_simple();
   bool IsConvex    = bot_poly.is_convex();
   if (IsConvex)
   {
   	std::cout<<"polygon is convex"<<std::endl;
   }
   if (!IsConvex)
   {
   	std::cout<<"polygon is nonconvex"<<std::endl;
   }
   // bool IsClockwise = (bot_poly.orientation() == CGAL::CLOCKWISE);
   // double Area      = bot_poly.area();
   // std::cout << "polygon p is";
   // if (!IsSimple) std::cout << " not";
   // std::cout << " simple." << std::endl;
   // std::cout << "polygon p is";
   // if (!IsConvex) std::cout << " not";
   // std::cout << " convex." << std::endl;
   // std::cout << "polygon p is";
   // if (!IsClockwise) std::cout << " not";
   // std::cout << " clockwise oriented." << std::endl;
   // std::cout << "the area of polygon p is " << Area << std::endl;
   // std::cout << std::endl;
if (!IsSimple)
   {
      std::cout<<"polygon is not simple go back to LABEL1"<<std::endl;
      delete [] disper_point;
      delete [] poly_seg;
      goto LABEL1;
      
   }
   // get edge line vector
   // for (int j = 0; j < EDGE_POLY; j++)
   // {
   //  if (j + 1 < EDGE_POLY)
   //   {edge_line[j].x = vertex[j].y-vertex[j+1].y;
   //        edge_line[j].y = vertex[j+1].x-vertex[j].x;
   //        edge_line[j].z = vertex[j].x *vertex[j+1].y - vertex[j+1].x*vertex[j].y;}
   //  if (j + 1 == EDGE_POLY)
   //  {
   //    edge_line[j].x = vertex[j].y-vertex[0].y;
   //    edge_line[j].y = vertex[0].x-vertex[j].x;
   //    edge_line[j].z = vertex[j].x *vertex[0].y - vertex[0].x*vertex[j].y;
   //  }

   // }
   // std::cout<<"get edge line done"<<std::endl;

   // sparse points 
   for (int j =0; j<50; j++)
   {
     for(int q =0; q<50; q++)
     {
      disper_point[50*j+q].x= -RADIUS+((2*RADIUS)/50)*q;
      disper_point[50*j+q].y= -RADIUS+((2*RADIUS)/50)*j;
      // std::cout<<"sparse points:"<<disper_point[50*j+q].x<<" "<<disper_point[50*j+q].y<<std::endl;
     }
     
   }

   // std::cout<<"sparse points done"<<std::endl;

   // get inside points
   vector <vector <double> > inside_points;
   for (int j = 0; j< 2500; j++)
   {
    bool IsInside = (bot_poly.bounded_side(Point_2(disper_point[j].x,disper_point[j].y)) == CGAL::ON_BOUNDED_SIDE);
    // std::cout<< "IsInside:"<<IsInside << std::endl;
    if (IsInside)
    {
      
      std::vector<double> one_inside_point;
      for(int k = 0; k < 2; k++)
      {
        if (k==0) one_inside_point.push_back(disper_point[j].x);
        if (k==1) one_inside_point.push_back(disper_point[j].y);
      }
      inside_points.push_back(one_inside_point);
      // std::cout<<"put back one_inside_point"<<std::endl;
    }
   }
   if (inside_points.size()==0) 
    {
      std::cout<<"no inside points go back to LABEL1"<<std::endl;
      // for ( int i = 0; i < EDGE_POLY; i++)
   	  // {
	      
	     //  std::cout<<"vertex"<<i<<":"<<"("<<vertex[i].x<<" ,"<<vertex[i].y<<" ,"<<vertex[i].z<<")"<<std::endl;
	     //  // std::cout<<"vertex"<<i<<":"<<"("<<vertex[i].x<<" ,"<<vertex[i].y<<" ,"<<vertex[i].z<<")"<<std::endl;
	     
      // }

      delete [] disper_point;
      delete [] poly_seg;
      goto LABEL1;
      
    }

   // std::cout<<"get inside points done"<<std::endl;
   // std::cout<<"get inside points vector size "<<inside_points.size()<<std::endl;
   for (int j= 0; j < inside_points.size(); j++)
   {
    bool IsContain = false;
    bool IsIntersect = false;
    Segment_2 grid_seg[4];
    Point_2 Plu(inside_points[j][0]-0.01,inside_points[j][1]+0.01);
    Point_2 Pld(inside_points[j][0]-0.01,inside_points[j][1]-0.01);
    Point_2 Pru(inside_points[j][0]+0.01,inside_points[j][1]+0.01);
    Point_2 Prd(inside_points[j][0]+0.01,inside_points[j][1]-0.01);
    grid_seg[0]= Segment_2(Plu,Pld);
    grid_seg[1]= Segment_2(Pld,Prd);
    grid_seg[2]= Segment_2(Prd,Pru);
    grid_seg[3]= Segment_2(Pru,Plu);

    IsContain = (bot_poly.bounded_side(Plu) == CGAL::ON_BOUNDED_SIDE) && (bot_poly.bounded_side(Pld) == CGAL::ON_BOUNDED_SIDE) && 
                (bot_poly.bounded_side(Pru) == CGAL::ON_BOUNDED_SIDE) && (bot_poly.bounded_side(Prd) == CGAL::ON_BOUNDED_SIDE);
    for (int k = 0; k<EDGE_POLY; k++ )
    {
      for (int e = 0; e < 4; e++)
      {
        CGAL::cpp11::result_of<Intersect_2(Segment_2, Segment_2)>::type
        result = intersection(grid_seg[e], poly_seg[k]);
        IsIntersect = IsIntersect || result;
        // if(IsIntersect)
        //   {
        //     delete [] disper_point;
        //     delete [] poly_seg;
        //     goto LABEL1;
        //   }
      }
    }
    if (IsContain && !IsIntersect) 
    	{
    		std::cout<<" grid inside go back to LABEL2"<<std::endl;
    		goto LABEL2;
    	}
    if ((!IsContain && j == inside_points.size()-1) || IsIntersect) 
    	{
			std::cout<<"no grid inside go back to LABEL1"<<std::endl;
			delete [] disper_point;
			delete [] poly_seg;
			goto LABEL1;
			

		}
   }
   // double** disper_dis = new double*[inside_points.size()]();
   // for(int j = 0; j < inside_points.size(); ++j) {
   //    disper_dis[j] = new (nothrow) double[EDGE_POLY]();
   // }
   // if (disper_dis == nullptr)
   // std::cout << "Error: disper_dis memory could not be allocated";

   // // get distance 
   // for (int j =0; j<inside_points.size(); j++)
   // {
   //   for(int k = 0; k < EDGE_POLY; k++)
   //   {
   //      disper_dis[j][k]= abs(edge_line[k].x * inside_points[j][0] + edge_line[k].y * inside_points[j][1] + edge_line[k].z) / sqrt(edge_line[k].x * edge_line[k].x + edge_line[k].y * edge_line[k].y);
   //   }
   // }

   // // std::cout<<"get inside points distance done"<<std::endl;
   // // std::cout<<"get inside points distance size "<<disper_dis[0][2]<<std::endl;

   // // get largest distance
   // double largest_incircle_radius = RADIUS;
   // for (int j = 0; j < inside_points.size(); j++)
   // {
   //  for (int k = 0; k < EDGE_POLY; k++)
   //  {
   //    bool distance_eql = true;
   //    for (int p = k + 1; p < EDGE_POLY; p++)
   //    {
   //      if (disper_dis[j][p]-disper_dis[j][k]>0.001) 
   //        {
   //          distance_eql= false;
   //          break;
   //        }
   //  }
   //  if (k== EDGE_POLY-1 && distance_eql == true) 
   //    for ( int q = 0; q< EDGE_POLY; q++)
   //    largest_incircle_radius = (disper_dis[j][q]<largest_incircle_radius)?disper_dis[j][k]:largest_incircle_radius;
   //  }
   // }
   // std::cout<<"get largest distance done"<<std::endl;
   // std::cout<<"get largest distance size "<<largest_incircle_radius<<std::endl;
   LABEL2:
   delete [] disper_point;
   delete [] poly_seg;
   // delete [] edge_line;
   // for(int j = 0; j < inside_points.size(); j++) {
   //   delete [] disper_dis[j];
   // }

   // if (largest_incircle_radius < 0.03 || bot_poly.area()< 0.02) goto LABEL1;
   // std::cout<<"get largest distance of poly"<<poly_no<< " "<<largest_incircle_radius<<std::endl;
   // double Area      = bot_poly.area();
   // if ( Area < 0.02)
   // std::cout << "the area of polygon"<< poly_no<<" is " << Area << std::endl;


   // std::cout<< "Scale Done"<<std::endl;
   // std::cout<< "Transfer vertex..."<<std::endl; 
   TRANS:
   // double transform_seed = ((double)rand() / double(RAND_MAX))*2+(-1);
   x_range_left = -bot_poly.bbox().xmin()-(RADIUS-0.01);
   x_range_right = RADIUS-0.01-bot_poly.bbox().xmax();
   y_range_top = RADIUS-0.01-bot_poly.bbox().ymax();
   y_range_bot = -bot_poly.bbox().ymin()-(RADIUS-0.01);
   delta_x = ((double)rand() / double(RAND_MAX))*(x_range_right-x_range_left)+(x_range_left);
   delta_y = ((double)rand() / double(RAND_MAX))*(y_range_top-y_range_bot)+(y_range_bot);
   Polygon_2 bot_poly_trans;
   bot_poly_trans.clear();
   // delta_x = x_range_left;
   // delta_y = y_range_top;
   for ( int i = 0; i < EDGE_POLY; i++)
   {
      vertex[i] = { vertex[i].x + delta_x, vertex[i].y + delta_y, 0 };
      std::cout<<"vertex"<<i<<":"<<"("<<vertex[i].x<<" ,"<<vertex[i].y<<" ,"<<vertex[i].z<<")"<<std::endl;
      // std::cout<<"vertex"<<i<<":"<<"("<<vertex[i].x<<" ,"<<vertex[i].y<<" ,"<<vertex[i].z<<")"<<std::endl;
      bot_poly_trans.push_back(Point_2(vertex[i].x, vertex[i].y ));
   }

   if (!bot_poly_trans.is_simple())
   {
    std::cout<<"bottom polygon is not simple go back to LABEL1"<<std::endl;
    goto LABEL1;
   }

   if(bot_poly_trans.bbox().xmax()>(RADIUS-0.01) || bot_poly_trans.bbox().ymax()>(RADIUS-0.01) || bot_poly_trans.bbox().xmin()<(-RADIUS+0.01) || bot_poly_trans.bbox().ymin()<(-RADIUS+0.01))
   {
    std::cout<<"Error, out of boundary:"<< bot_poly_trans.bbox()<<endl;
    std::cout<<"vertex number: "<<bot_poly_trans.size()<<endl;
    goto TRANS;
   }
   std::cout<<"trans done"<<std::endl;


   // std::cout<< "Transfer Done"<<std::endl;
   
   for ( int  i = EDGE_POLY; i< 2*EDGE_POLY; i++) {vertex[i] = vertex[i-EDGE_POLY]+ extrude;}
      
   // define top-bottom face
   int p = 0;
   for ( int i = 0; i < EDGE_POLY; i++) {top_bott[i][p]=i;}
   p = 1;
   for ( int i = 0; i < EDGE_POLY; i++) {i == 0 ? top_bott[i][p]=EDGE_POLY : top_bott[i][p]=2*EDGE_POLY-i;}   

   // define side face
   int i = 0;
   for ( p = 0; p < EDGE_POLY; p++ ) {face[p][i] = p;}
   i = 1;
   for ( p = 0; p < EDGE_POLY; p++ ) {face[p][i] = face[p][0] + EDGE_POLY;}
   i = 2;
   for ( p = 0; p < EDGE_POLY; p++ ) {EDGE_POLY+1+p > 2*EDGE_POLY-1 ? face[p][i] = EDGE_POLY : face[p][i] = EDGE_POLY+1+p ;}
   i = 3;
   for ( p = 0; p < EDGE_POLY; p++ ) {1 + p > EDGE_POLY-1 ? face[p][i] = 0 : face[p][i] = 1 + p ;}

   
   

}

void output( int mesh_no, int EDGE_POLY, Pt* vertex, int** face, int** top_bott, Surface_mesh &mesh )
{
   #define SF <<fixed<< setprecision( 6 ) <<
   #define SI <<setw( 5 ) <<
   ofstream outfiles;
   string filename;
   stringstream a;
   a << mesh_no;
   filename = "mesh" + a.str() + ".off";
   outfiles.open(filename.c_str(), ios::out);
   if (outfiles.is_open())
   {
      outfiles << "OFF" << endl;
      outfiles << 2*EDGE_POLY SI EDGE_POLY+2 SI "0"<< endl;
      // cout << "VERTICES:\n";
      for ( int v = 0; v < 2 * EDGE_POLY; v++ ) outfiles SF vertex[v].x << "  "  SF vertex[v].y << "  "  SF vertex[v].z << '\n';

      // cout << "\nFACE:\n";
      // output top and bottom face 
      for ( int p = 0; p < 2; p++ ) 
      {
         outfiles << EDGE_POLY;
         for ( int i = 0; i < EDGE_POLY; i++ ) outfiles SI top_bott[i][p] << " ";
         outfiles << '\n';
      }
      // output side face
      for ( int p = 0; p < EDGE_POLY; p++ ) 
      {
         outfiles << "4";
         for ( int i = 0; i < 4; i++ ) outfiles SI face[p][i] << " ";
         outfiles << '\n';
      }
   }
   outfiles.close();

   std::ifstream input(filename);
    // Surface_mesh mesh;

    if (!input || !(input >> mesh) || mesh.is_empty())//
    {
      std::cerr << "Not a valid off file/output." << std::endl;
      
    }
    input.close();

    CGAL::Polygon_mesh_processing::triangulate_faces(mesh);

    std::remove(filename.c_str());
}

void off2obj(int mesh_no,string output_file_path)   //convert off file to obj file
{

  string infiles;
  string outfiles;
  stringstream a;
  a << mesh_no;
  infiles = "merged_mesh.off";
  istream*     p_in  = &cin;
  ifstream     in;
  
  in.open( infiles);
  p_in = &in;
  
  
  if ( !*p_in) {
      cerr  << ": error: cannot open file '"<< infiles
           << "' for reading." <<endl;
      
  }

  outfiles = output_file_path + "poly_"+ a.str()+".obj";
  ostream*     p_out = &cout;
  ofstream     out;
  
  out.open( outfiles);
  p_out = &out;
  
  
  if ( !*p_out) {
      cerr  << ": error: cannot open file '"<< outfiles
           << "' for writing." <<endl;
      
  }

  
  CGAL::File_writer_wavefront  writer;
  CGAL::generic_copy_OFF( *p_in, *p_out, writer);
  

  if ( !*p_in) {
      cerr  << " read error: while reading file '"<< infiles << "'."
           << endl;
      
  }
  if ( !*p_out) {
      cerr  <<" write error: while writing file '"<< outfiles << "'."
           << endl;
      
  }
  
  std::remove(infiles.c_str());
}

void create_ground ( double length, double width, double height, double target_edge_length, Surface_mesh &mesh)
{
	std::ifstream input("cube.off");
    

    if (!input || !(input >> mesh) || mesh.is_empty())//
    {
      std::cerr << "Not a valid off file/iso." << std::endl;
      
    }
    input.close();
    double scale_x = length, scale_y = width, scale_z = height;
    Aff_transformation_3 transform(RT_num(scale_x), RT_num(0),      RT_num(0),      RT_num(0),
                                   RT_num(0),       RT_num(scale_y),RT_num(0),      RT_num(0),
                                   RT_num(0),       RT_num(0),      RT_num(scale_z),RT_num(-0.005-0.5*scale_z),
                                                                                    RT_num(1));
    PMP::transform (transform, mesh); 
    //CGAL::Polygon_mesh_processing::triangulate_faces(mesh);
    ISO_MESH(1, target_edge_length, mesh,999);

}