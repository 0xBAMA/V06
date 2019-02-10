#include "../voraldo/v.h"

using std::cout;
using std::endl;

//---------------------------
Voraldo_IO::Voraldo_IO(Voraldo *p)
{
 parent = p;
}

Voraldo_IO::~Voraldo_IO()
{

}

void Voraldo_IO::display(std::string filename, double x_rot, double y_rot, double z_rot, double scale, bool perspective)
{
  using namespace cimg_library;

 	int image_x_dimension = 1921;
 	int image_y_dimension = 1081;

 	CImg<unsigned char> img(image_x_dimension,image_y_dimension,1,3,0);

 	const unsigned char	gold[3] = {255,215,0};
 	const unsigned char dark_gold[3] = {127,107,0};
 	const unsigned char white[3] = {255,255,255};
 	const unsigned char black[3] = {0,0,0};
 	const unsigned char pink[3] = {255,0,255};


 	//frame top
 	img.draw_line(0,1,image_x_dimension,1,gold);
 	img.draw_line(0,2,image_x_dimension,2,dark_gold);

 	//frame bottom
 	img.draw_line(0,image_y_dimension-2,image_x_dimension,image_y_dimension-2,dark_gold);
 	img.draw_line(0,image_y_dimension-1,image_x_dimension,image_y_dimension-1,gold);

 	//frame left from low x low y to low x high y
 	img.draw_line(1,0,1,image_y_dimension,gold);
 	img.draw_line(2,0,2,image_y_dimension,dark_gold);

 	//frame right from high x low y to high x high y
 	img.draw_line(image_x_dimension-2,0,image_x_dimension-2,image_y_dimension,dark_gold);
 	img.draw_line(image_x_dimension-2,0,image_x_dimension-2,image_y_dimension,gold);


 	int block_xdim = parent->x_dim;
 	int block_ydim = parent->y_dim;
 	int block_zdim = parent->z_dim;

 	int block_xdim_squared = block_xdim*block_xdim;
 	int block_ydim_squared = block_ydim*block_ydim;
 	int block_zdim_squared = block_zdim*block_zdim;

 	vec d_center = vec(block_xdim/2.0,block_ydim/2.0,block_zdim/2.0);

 	mat rotation_x_axis;
 //refernces [column][row]
 		rotation_x_axis[0][0] = 1;
 		rotation_x_axis[0][1] = 0;
 		rotation_x_axis[0][2] = 0;
 		rotation_x_axis[1][0] = 0;
 		rotation_x_axis[1][1] = std::cos(x_rot);
 		rotation_x_axis[1][2] = std::sin(x_rot);
 		rotation_x_axis[2][0] = 0;
 		rotation_x_axis[2][1] = -1.0*std::sin(x_rot);
 		rotation_x_axis[2][2] = std::cos(x_rot);

 	mat rotation_y_axis;
 		rotation_y_axis[0][0] = std::cos(y_rot);
 		rotation_y_axis[0][1] = 0;
 		rotation_y_axis[0][2] = -1.0*std::sin(y_rot);
 		rotation_y_axis[1][0] = 0;
 		rotation_y_axis[1][1] = 1;
 		rotation_y_axis[1][2] = 0;
 		rotation_y_axis[2][0] = std::sin(y_rot);
 		rotation_y_axis[2][1] = 0;
 		rotation_y_axis[2][2] = std::cos(y_rot);

 	mat rotation_z_axis;
 		rotation_z_axis[0][0] = std::cos(z_rot);
 		rotation_z_axis[0][1] = std::sin(z_rot);
 		rotation_z_axis[0][2] = 0;
 		rotation_z_axis[1][0] = -1.0*std::sin(z_rot);
 		rotation_z_axis[1][1] = std::cos(z_rot);
 		rotation_z_axis[1][2] = 0;
 		rotation_z_axis[2][0] = 0;
 		rotation_z_axis[2][1] = 0;
 		rotation_z_axis[2][2] = 1;

 	vec d_xvec = mul(rotation_x_axis,mul(rotation_y_axis,mul(rotation_z_axis,vec(1,0,0))));
 	vec d_yvec = mul(rotation_x_axis,mul(rotation_y_axis,mul(rotation_z_axis,vec(0,1,0))));
 	vec d_zvec = mul(rotation_x_axis,mul(rotation_y_axis,mul(rotation_z_axis,vec(0,0,1))));

 	//tip radius is the radius of the sphere that touches the tips of the block's corners
 	double tip_radius = std::sqrt(block_xdim_squared+block_ydim_squared+block_zdim_squared)/2.0;
 	double max_dist = 2*2.2*tip_radius;
 	double min_dist = 0.2*tip_radius;
 	//factor of two is for the incremental step of length 0.5
 	//factor of 2.2 is for the tip radius plus the camera sphere radius, 1+1.2

 	vec cam_position = d_center - 1.2*tip_radius*d_yvec;

 	vec cam_up = scale*d_zvec; //may need to change the scaling
 	vec cam_right = scale*d_xvec;



 	int image_center_x = (image_x_dimension-1)/2;;
 	int image_center_y = (image_y_dimension-1)/2;;

 	int image_current_x, image_current_y;

 	vec vector_starting_point;
 	vec vector_increment = 0.5*normalize(-1.0*(cam_position-d_center));

 	vec vector_test_point;

 	vec block_min = vec(0,0,0);
 	vec block_max = vec(block_xdim,block_ydim,block_zdim);

 	Vox tempstate;
 	unsigned char	point_color[3];

 	bool xtest,ytest,ztest;
 	bool line_box_intersection;

 	//double t0 = 0;
 	//double t1 = 9999;

  double tmin, tmax;

 	double tintersect; //for line/box intersection

 	for(double x = -(image_x_dimension/2-5); x <= (image_x_dimension/2-5); x++)
 		for(double y = -(image_y_dimension/2-5); y <= (image_y_dimension/2-5); y++)
 		{//init (reset)
 			line_box_intersection = false;

 			image_current_x = image_center_x + x;
 			image_current_y = image_center_y + y;

 			// if(perspective == true) //this gets added inside the loop - note that the linetest will have to consider the perspective corrected ray
 			// 	vector_increment_perspective = vector_increment + x*0.1*cam_right - y*0.1*cam_up;

 			vector_starting_point = cam_position + x*cam_up + y*cam_right;

 			//figure out if the parametric line established by parameter z and
 			//	L = vector_starting_point + z*vector_increment
 			//intersects with the box established by (0,0,0) (x,y,z)
 			//i.e. block_min and block_max

 			//The goal here is to achieve some level of speedup when compared to
 			//doing an exhaustive check through all points on all vectors from
 			//the pixels in the image.
 			line_box_intersection = parent->intersect_ray_bbox(block_min,block_max,vector_starting_point,vector_increment,tmin,tmax);

 			if(line_box_intersection)
 			{//we will need to step through the box
 				for(double z = tmin; z <= tmax; z+=0.5) //go from close intersection point to far intersection point
 				{
 					vector_test_point = vector_starting_point + z*vector_increment;
 					tempstate = parent->get_data_by_vector_index(vector_test_point);
 					if(tempstate.state!=0)
 					{
 						point_color = palette[temp];

 						img.draw_point(image_current_x,image_current_y,point_color);
 						break;
 					}
 				}
 			}
 			else
 			{
 				img.draw_point(image_current_x,image_current_y,black);
 			}
 		}

 	img.save_bmp(filename.c_str());
}




//---------------------------
Voraldo_Draw::Voraldo_Draw(Voraldo *p)
{
 parent = p;
}

Voraldo_Draw::~Voraldo_Draw()
{

}

void Voraldo_Draw::init_block(int x, int y, int z)
{

 if(parent->data != NULL)
 {
		delete[] parent->data;
 }

 parent->x_dim = x;
	parent->y_dim = y;
	parent->z_dim = z;

	parent->num_cells = parent->x_dim * parent->y_dim * parent->z_dim;

 parent->data = new Vox[parent->num_cells];

 for(int i = 0; i < parent->num_cells; i++)
 {
  parent->data[i].mask = false;
  parent->data[i].state = 0;
  parent->data[i].alpha = 0;
 }
}

void Voraldo_Draw::clear_all()
{
   for(int i = 0; i < parent->num_cells; i++)
   {
     if(!parent->data[i].mask)
     {
      parent->data[i].state = 0;
      parent->data[i].alpha = 0;
      //deliberately do not unmask - separate function for that operation
     }
   }
}

void Voraldo_Draw::mask_unmask_all()
{
 for(int i = 0; i < parent->num_cells; i++)
 {
  parent->data[i].mask = false;
 }
}

void Voraldo_Draw::mask_invert_mask()
{
 for(int i = 0; i < parent->num_cells; i++)
 {
  parent->data[i].mask = parent->data[i].mask ? false:true;
 }
}

void Voraldo_Draw::mask_all_nonzero()
{
 for(int i = 0; i < parent->num_cells; i++)
 {
  if(parent->data[i].state != 0)
  {
   parent->data[i].mask = true;
  }
 }
}

void Voraldo_Draw::mask_by_state(unsigned char s)
{
 for(int i = 0; i < parent->num_cells; i++)
 {
  if(parent->data[i].state == s)
  {
   parent->data[i].mask = true;
  }
 }
}

void Voraldo_Draw::draw_noise(bool draw, unsigned char alpha, bool mask)
{
  for(int i = 0; i < parent->num_cells; i++)
  {
     if(std::rand()%696 == 69)
     {
       if(!parent->data[i].mask)
       {
        parent->data[i].state = std::rand()%255;//this is a little different
        parent->data[i].alpha = alpha;
        parent->data[i].mask = mask;
       }
     }
  }
}

void Voraldo_Draw::draw_point(vec point, Vox set, bool draw, bool mask)
{
 parent->set_data_by_vector_index(point,set,draw,mask);
}

void Voraldo_Draw::draw_line_segment(vec v1, vec v2, Vox set, bool draw, bool mask)
{
 vec starting_point = v1;
	vec current_point = starting_point;
	vec line_vector = (v2-v1);

	int length = std::floor(linalg::length(line_vector));

	for(int i = 0; i < length; i++)
	{
		current_point[0] = starting_point[0] + i*(line_vector[0]/length);
		current_point[1] = starting_point[1] + i*(line_vector[1]/length);
		current_point[2] = starting_point[2] + i*(line_vector[2]/length);
		draw_point(current_point,set);
	}
}

void Voraldo_Draw::draw_triangle(vec v0, vec v1, vec v2, Vox set, bool draw, bool mask)
{
 //point zero is the origin point

 	vec side1 = v1-v0;
 	vec side2 = v2-v0;

 	vec c1(0,0,0);
 	vec c2(0,0,0);
  //c1 and c2 are the ends of the line segment that is going to be drawn
  // these endpoints are ranged along the two sides so as to achieve coverage

 	double length;

 	if(linalg::length(side1) > linalg::length(side2))
 	{
 		length = std::floor(linalg::length(side1));
 	}
 	else
 	{
 		length = std::floor(linalg::length(side2));
 	}

 	if(length <= 2)
  {
 		draw_point(v0,set);
 		draw_point(v1,set);
 		draw_point(v2,set);
 	}
  else
  {

 		side1 = side1/length;
 		side2 = side2/length;

 		for(int i = 0; i < length; i++)
 		{
 			c1[0] = v0[0] + i*side1[0];
 			c1[1] = v0[1] + i*side1[1];
 			c1[2] = v0[2] + i*side1[2];

 			c2[0] = v0[0] + i*side2[0];
 			c2[1] = v0[1] + i*side2[1];
 			c2[2] = v0[2] + i*side2[2];

 			draw_line_segment(c1,c2,set);
 		}
 	}
}

void Voraldo_Draw::draw_sphere(vec center_point, double radius, Vox set, bool draw, bool mask)
{
 vec index;
 for(int i = 0; i < parent->x_dim; i++)
 	{
 		for(int j = 0; j < parent->y_dim; j++)
 		{
 			for(int k = 0; k < parent->z_dim; k++)
 			{
 				double testx = (i-center_point[0])*(i-center_point[0]);	//apply offsets and square
 				double testy = (j-center_point[1])*(j-center_point[1]);
 				double testz = (k-center_point[2])*(k-center_point[2]);

 				if((testx + testy + testz) < radius*radius)
 				{	//pretty self explainatory, equation of sphere
      index = vec(i,j,k);
 					parent->set_data_by_vector_index(index,set,draw,mask);
 				}
 			}
 		}
 	}
}

void Voraldo_Draw::draw_ellipsoid(vec center_point, vec radii, Vox set, bool draw, bool mask)
{
 vec index;
 for(int i = 0; i < parent->x_dim; i++)
 	{
 		for(int j = 0; j < parent->y_dim; j++)
 		{
 			for(int k = 0; k < parent->z_dim; k++)
 			{
 				double testx = (i-center_point[0])*(i-center_point[0]);	//apply offsets and square
 				double testy = (j-center_point[1])*(j-center_point[1]);
 				double testz = (k-center_point[2])*(k-center_point[2]);

 				double radx = radii[0]*radii[0];
 				double rady = radii[1]*radii[1];
 				double radz = radii[2]*radii[2];

 				double result = testx/radx + testy/rady + testz/radz;

 				if(result <= 1){	//point inside ellipsoid - do we want to be able to invert this?
 					//(outside, or on the surface, with >= and ==, respectively)
      index = vec(i,j,k);
 					parent->set_data_by_vector_index(index,set,draw,mask);
 				}
 			}
 		}
 	}
}

void Voraldo_Draw::draw_cylinder(vec bvec, vec tvec, double radius, Vox set, bool draw, bool mask)
{
 vec ndirection = tvec - bvec;

	auto bx0 = bvec[0]; auto ba = ndirection[0];	auto tx0 = tvec[0]; auto ta = ndirection[0];
	auto by0 = bvec[1]; auto bb = ndirection[1];	auto ty0 = tvec[1]; auto tb = ndirection[1];
	auto bz0 = bvec[2]; auto bc = ndirection[2];	auto tz0 = tvec[2]; auto tc = ndirection[2];

	//I did this on a whiteboard

	double bplanetest = 0.0;
	double tplanetest = 0.0;

	double point_to_line_distance = 0.0;

 vec index;

	for(int i = 0; i < parent->x_dim; i++){
		for(int j = 0; j < parent->y_dim; j++){
			for(int k = 0; k < parent->z_dim; k++){
				//planetests
				bplanetest = ba * (i - bx0) + bb * (j - by0) + bc * (k - bz0);
				tplanetest = ta * (i - tx0) + tb * (j - ty0) + tc * (k - tz0);

				//using the basic equation for a plane, we can do an interesting test

				//These variables will be greater than zero if the test point is on the side of the plane
				//that the normal vector is pointing towards, and less than zero if the test point is on
				//the side of the plane that the normal vector is not pointing towards. That is to say, in
				//my case - bplanetest tells me whether the point is above the bottom plane, and tplanetest
				//tells me whether the point is above the top plane. If it is above the bottom, and below
				//the top - we are within the ends of the cylinder, and can proceed. Thus, the condition
				//for the following if statement:

				if(bplanetest >= 0 && tplanetest <= 0){
					//do the point to line distance thing
					//algorithm from http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html

					//This test takes as a precondition that the point is determined to be within the slice
					//of space established by two parallel planes. Now that this is known, we will see how
					//far they are from the line, which runs perpendicular to both planes. If both of these
					//sequential tests evaluate positively, we know we are within the extents of the cylinder.

					point_to_line_distance = linalg::length(cross(tvec-bvec,bvec-vec(i,j,k)))/linalg::length(tvec-bvec);
					if(point_to_line_distance <= radius){
						parent->set_data_by_vector_index(index,set,draw,mask);
					}
				}
			}
		}
	}
}

void Voraldo_Draw::draw_tube(vec bvec, vec tvec, double inner_radius, double outer_radius, Vox set, bool draw, bool mask)
{
 vec ndirection = tvec - bvec;

 	auto bx0 = bvec[0]; auto ba = ndirection[0];	auto tx0 = tvec[0]; auto ta = ndirection[0];
 	auto by0 = bvec[1]; auto bb = ndirection[1];	auto ty0 = tvec[1]; auto tb = ndirection[1];
 	auto bz0 = bvec[2]; auto bc = ndirection[2];	auto tz0 = tvec[2]; auto tc = ndirection[2];

 	//I did this on a whiteboard

 	double bplanetest = 0.0;
 	double tplanetest = 0.0;

 	double point_to_line_distance = 0.0;
  vec index;

 	for(int i = 0; i < parent->x_dim; i++){
 		for(int j = 0; j < parent->y_dim; j++){
 			for(int k = 0; k < parent->z_dim; k++){
 				//planetests
 				bplanetest = ba * (i - bx0) + bb * (j - by0) + bc * (k - bz0);
 				tplanetest = ta * (i - tx0) + tb * (j - ty0) + tc * (k - tz0);

 				//using the basic equation for a plane, we can do an interesting test

 				//These variables will be greater than zero if the test point is on the side of the plane
 				//that the normal vector is pointing towards, and less than zero if the test point is on
 				//the side of the plane that the normal vector is not pointing towards. That is to say, in
 				//my case - bplanetest tells me whether the point is above the bottom plane, and tplanetest
 				//tells me whether the point is above the top plane. If it is above the bottom, and below
 				//the top - we are within the ends of the cylinder, and can proceed. Thus, the condition
 				//for the following if statement:

 				if(bplanetest >= 0 && tplanetest <= 0){
 					//do the point to line distance thing
 					//algorithm from http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html

 					//This test takes as a precondition that the point is determined to be within the slice
 					//of space established by two parallel planes. Now that this is known, we will see how
 					//far they are from the line, which runs perpendicular to both planes. If both of these
 					//sequential tests evaluate positively, we know we are within the extents of the cylinder.

 					point_to_line_distance = linalg::length(cross(tvec-bvec,bvec-vec(i,j,k)))/linalg::length(tvec-bvec);

 					if(point_to_line_distance <= outer_radius && point_to_line_distance >= inner_radius){
       parent->set_data_by_vector_index(index,set,draw,mask);
 					}
 				}
 			}
 		}
 	}

}

void Voraldo_Draw::draw_quadrilateral_hexahedron(vec a, vec b, vec c, vec d, vec e, vec f, vec g, vec h, Vox set, bool draw, bool mask)
{
 vec center = a + b + c + d + e + f + g + h;
 	center = vec(center[0]/8, center[1]/8, center[2]/8);

 	bool plusx1 = false;
 	bool plusx2 = false;

 //    888   Y88b  /
 //  __888__  Y88b/
 //    888     Y88b
 //    888     /Y88b
 //           /  Y88b

 	//CDGH

 	//TRIANGLES ARE CDH, CGH

 	vec plusx_side1 = vec(c[0]-d[0],c[1]-d[1],c[2]-d[2]);
 	vec plusx_side2 = vec(h[0]-d[0],h[1]-d[1],h[2]-d[2]);

 	vec plusx_side3 = vec(c[0]-g[0],c[1]-g[1],c[2]-g[2]);
 	vec plusx_side4 = vec(h[0]-g[0],h[1]-g[1],h[2]-g[2]);

 	vec plusx_planevec1 = cross(plusx_side1,plusx_side2);
 	vec plusx_planept1 = d;

 	vec plusx_planevec2 = cross(plusx_side3,plusx_side4);
 	vec plusx_planept2 = g;

 	if(!parent->planetest(plusx_planept1,plusx_planevec1,center))
 	{	//if the center fails the plane test, flip the normal
 		plusx_planevec1 = vec(-1*plusx_planevec1[0],-1*plusx_planevec1[1],-1*plusx_planevec1[2]);
 	}

 	if(!parent->planetest(plusx_planept2,plusx_planevec2,center))
 	{	//if the center fails the plane test, flip the normal
 		plusx_planevec2 = vec(-1*plusx_planevec2[0],-1*plusx_planevec2[1],-1*plusx_planevec2[2]);
 	}

 //-------------------------------------------


 	bool minusx1 = false;
 	bool minusx2 = false;

 //       Y88b  /
 //  ____  Y88b/
 //         Y88b
 //         /Y88b
 //        /  Y88b

 	//ABEF

 	//TRIANGLES ARE ABF, AEF

 	vec minusx_side1 = vec(a[0]-b[0],a[1]-b[1],a[2]-b[2]);
 	vec minusx_side2 = vec(f[0]-b[0],f[1]-b[1],f[2]-b[2]);

 	vec minusx_side3 = vec(a[0]-e[0],a[1]-e[1],a[2]-e[2]);
 	vec minusx_side4 = vec(f[0]-e[0],f[1]-e[1],f[2]-e[2]);

 	vec minusx_planevec1 = cross(minusx_side1,minusx_side2);
 	vec minusx_planept1 = b;

 	vec minusx_planevec2 = cross(minusx_side3,minusx_side4);
 	vec minusx_planept2 = e;

 	if(!parent->planetest(minusx_planept1,minusx_planevec1,center))
 	{	//if the center fails the plane test, flip the normal
 		minusx_planevec1 = vec(-1*minusx_planevec1[0],-1*minusx_planevec1[1],-1*minusx_planevec1[2]);
 	}

 	if(!parent->planetest(minusx_planept2,minusx_planevec2,center))
 	{	//if the center fails the plane test, flip the normal
 		minusx_planevec2 = vec(-1*minusx_planevec2[0],-1*minusx_planevec2[1],-1*minusx_planevec2[2]);
 	}

 //-------------------------------------------


 	bool plusy1 = false;
 	bool plusy2 = false;

 //    888   Y88b  /
 //  __888__  Y888/
 //    888     Y8/
 //    888      Y
 //            /
 //          _/

 	//ACEG

 	//TRIANGLES ARE ACG, AEG

 	vec plusy_side1 = vec(a[0]-c[0],a[1]-c[1],a[2]-c[2]);
 	vec plusy_side2 = vec(g[0]-c[0],g[1]-c[1],g[2]-c[2]);

 	vec plusy_side3 = vec(a[0]-e[0],a[1]-e[1],a[2]-e[2]);
 	vec plusy_side4 = vec(g[0]-e[0],g[1]-e[1],g[2]-e[2]);

 	vec plusy_planevec1 = cross(plusy_side1,plusy_side2);
 	vec plusy_planept1 = c;

 	vec plusy_planevec2 = cross(plusy_side3,plusy_side4);
 	vec plusy_planept2 = e;

 	if(!parent->planetest(plusy_planept1,plusy_planevec1,center))
 	{	//if the center fails the plane test, flip the normal
 		plusy_planevec1 = vec(-1*plusy_planevec1[0],-1*plusy_planevec1[1],-1*plusy_planevec1[2]);
 	}

 	if(!parent->planetest(plusy_planept2,plusy_planevec2,center))
 	{	//if the center fails the plane test, flip the normal
 		plusy_planevec2 = vec(-1*plusy_planevec2[0],-1*plusy_planevec2[1],-1*plusy_planevec2[2]);
 	}

 //-------------------------------------------


 	bool minusy1 = false;
 	bool minusy2 = false;

 //       Y88b  /
 //  ____  Y888/
 //         Y8/
 //          Y
 //         /
 //       _/

 	//BDFH

 	//TRIANGLES ARE BDH, BFH
 	vec minusy_side1 = vec(b[0]-d[0],b[1]-d[1],b[2]-d[2]);
 	vec minusy_side2 = vec(h[0]-d[0],h[1]-d[1],h[2]-d[2]);

 	vec minusy_side3 = vec(b[0]-f[0],b[1]-f[1],b[2]-f[2]);
 	vec minusy_side4 = vec(h[0]-f[0],h[1]-f[1],h[2]-f[2]);

 	vec minusy_planevec1 = cross(minusy_side1,minusy_side2);
 	vec minusy_planept1 = d;

 	vec minusy_planevec2 = cross(minusy_side3,minusy_side4);
 	vec minusy_planept2 = f;

 	if(!parent->planetest(minusy_planept1,minusy_planevec1,center))
 	{	//if the center fails the plane test, flip the normal
 		minusy_planevec1 = vec(-1*minusy_planevec1[0],-1*minusy_planevec1[1],-1*minusy_planevec1[2]);
 	}

 	if(!parent->planetest(minusy_planept2,minusy_planevec2,center))
 	{	//if the center fails the plane test, flip the normal
 		minusy_planevec2 = vec(-1*minusy_planevec2[0],-1*minusy_planevec2[1],-1*minusy_planevec2[2]);
 	}


 //-------------------------------------------


 	bool plusz1 = false;
 	bool plusz2 = false;


 //    888    ~~~d88P
 //  __888__    d88P
 //    888     d88P
 //    888    d88P
 //          d88P___

 	//ABCD

 	//TRIANGLES ARE ABD, ACD

 	vec plusz_side1 = vec(a[0]-b[0],a[1]-b[1],a[2]-b[2]);
 	vec plusz_side2 = vec(d[0]-b[0],d[1]-b[1],d[2]-b[2]);

 	vec plusz_side3 = vec(a[0]-c[0],a[1]-c[1],a[2]-c[2]);
 	vec plusz_side4 = vec(d[0]-c[0],d[1]-c[1],d[2]-c[2]);

 	vec plusz_planevec1 = cross(plusz_side1,plusz_side2);
 	vec plusz_planept1 = b;

 	vec plusz_planevec2 = cross(plusz_side3,plusz_side4);
 	vec plusz_planept2 = c;

 	if(!parent->planetest(plusz_planept1,plusz_planevec1,center))
 	{	//if the center fails the plane test, flip the normal
 		plusz_planevec1 = vec(-1*plusz_planevec1[0],-1*plusz_planevec1[1],-1*plusz_planevec1[2]);
 	}

 	if(!parent->planetest(plusz_planept2,plusz_planevec2,center))
 	{	//if the center fails the plane test, flip the normal
 		plusz_planevec2 = vec(-1*plusz_planevec2[0],-1*plusz_planevec2[1],-1*plusz_planevec2[2]);
 	}

 //-------------------------------------------


 	bool minusz1 = false;
 	bool minusz2 = false;

 //        ~~~d88P
 //  ____    d88P
 //         d88P
 //        d88P
 //       d88P___

 	//EFGH

 	//TRIANGLES ARE EFH, EGH
 	vec minusz_side1 = vec(e[0]-f[0],e[1]-f[1],e[2]-f[2]);
 	vec minusz_side2 = vec(h[0]-f[0],h[1]-f[1],h[2]-f[2]);

 	vec minusz_side3 = vec(e[0]-g[0],e[1]-g[1],e[2]-g[2]);
 	vec minusz_side4 = vec(h[0]-g[0],h[1]-g[1],h[2]-g[2]);

 	vec minusz_planevec1 = cross(minusz_side1,minusz_side2);
 	vec minusz_planept1 = f;

 	vec minusz_planevec2 = cross(minusz_side3,minusz_side4);
 	vec minusz_planept2 = g;

 	if(!parent->planetest(minusz_planept1,minusz_planevec1,center))
 	{	//if the center fails the plane test, flip the normal
 		minusz_planevec1 = vec(-1*minusz_planevec1[0],-1*minusz_planevec1[1],-1*minusz_planevec1[2]);
 	}

 	if(!parent->planetest(minusz_planept2,minusz_planevec2,center))
 	{	//if the center fails the plane test, flip the normal
 		minusz_planevec2 = vec(-1*minusz_planevec2[0],-1*minusz_planevec2[1],-1*minusz_planevec2[2]);
 	}


 //-------------------------------------------

 //  ╔╦╗┌─┐┌─┐┌┬┐  ╦  ┌─┐┌─┐┌─┐
 //   ║ ├┤ └─┐ │   ║  │ ││ │├─┘
 //   ╩ └─┘└─┘ ┴   ╩═╝└─┘└─┘┴

 	vec current;

 	for(int i = 0; i < parent->x_dim; i++)
 	{
 		for(int j = 0; j < parent->y_dim; j++)
 		{
 			for(int k = 0; k < parent->z_dim; k++)
 			{

 				current = vec(i,j,k);

 				bool plusxtest = parent->planetest(plusx_planept1,plusx_planevec1,current)&&parent->planetest(plusx_planept2,plusx_planevec2,current);
 				bool minusxtest = parent->planetest(minusx_planept1,minusx_planevec1,current)&&parent->planetest(minusx_planept2,minusx_planevec2,current);
 				bool plusytest = parent->planetest(plusy_planept1,plusy_planevec1,current)&&parent->planetest(plusy_planept2,plusy_planevec2,current);
 				bool minusytest = parent->planetest(minusy_planept1,minusy_planevec1,current)&&parent->planetest(minusy_planept2,minusy_planevec2,current);
 				bool plusztest = parent->planetest(plusz_planept1,plusz_planevec1,current)&&parent->planetest(plusz_planept2,plusz_planevec2,current);
 				bool minusztest = parent->planetest(minusz_planept1,minusz_planevec1,current)&&parent->planetest(minusz_planept2,minusz_planevec2,current);

 				bool xtest = plusxtest&&minusxtest;
 				bool ytest = plusytest&&minusytest;
 				bool ztest = plusztest&&minusztest;


 				if(xtest && ytest && ztest)
 				{
 					parent->set_data_by_vector_index(current,set,draw,mask);
 				}
 			}
 		}
 	}
}
