#include "../voraldo/v.h"

using std::cout;
using std::endl;


Vox get_vox(unsigned char state, unsigned char alpha, bool mask)
{
  Vox temp;

  temp.state = state;
  temp.alpha = alpha;
  temp.mask = mask;

  return temp;
}


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

 	int image_x_dimension = 1919;
 	int image_y_dimension = 1079;

 	CImg<unsigned char> img(image_x_dimension,image_y_dimension,1,3,0);

 	const unsigned char	gold[3] = {255,215,0};
 	const unsigned char dark_gold[3] = {127,107,0};
 	const unsigned char white[3] = {255,255,255};
 	const unsigned char black[3] = {0,0,0};
 	const unsigned char pink[3] = {255,0,255};

//draw the border for the image
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


 	int block_xdim = parent->x_dim;                 int block_ydim = parent->y_dim;                 int block_zdim = parent->z_dim;
 	int block_xdim_squared = block_xdim*block_xdim; int block_ydim_squared = block_ydim*block_ydim;	int block_zdim_squared = block_zdim*block_zdim;

 	vec d_center = vec(block_xdim/2.0,block_ydim/2.0,block_zdim/2.0);

//rotation matricies allowing rotation of the viewer's position
 	mat rotation_x_axis;
 //refernces [column][row] - sin and cos take arguments in radians
 		rotation_x_axis[0][0] = 1;                       rotation_x_axis[1][0] = 0;                      rotation_x_axis[2][0] = 0;
 		rotation_x_axis[0][1] = 0;                       rotation_x_axis[1][1] = std::cos(x_rot);        rotation_x_axis[2][1] = -1.0*std::sin(x_rot);
 		rotation_x_axis[0][2] = 0;                       rotation_x_axis[1][2] = std::sin(x_rot);        rotation_x_axis[2][2] = std::cos(x_rot);

 	mat rotation_y_axis;
 		rotation_y_axis[0][0] = std::cos(y_rot);         rotation_y_axis[1][0] = 0;                      rotation_y_axis[2][0] = std::sin(y_rot);
 		rotation_y_axis[0][1] = 0;                       rotation_y_axis[1][1] = 1;                      rotation_y_axis[2][1] = 0;
 		rotation_y_axis[0][2] = -1.0*std::sin(y_rot);    rotation_y_axis[1][2] = 0;                      rotation_y_axis[2][2] = std::cos(y_rot);

 	mat rotation_z_axis;
 		rotation_z_axis[0][0] = std::cos(z_rot);         rotation_z_axis[1][0] = -1.0*std::sin(z_rot);   rotation_z_axis[2][0] = 0;
 		rotation_z_axis[0][1] = std::sin(z_rot);         rotation_z_axis[1][1] = std::cos(z_rot);        rotation_z_axis[2][1] = 0;
 		rotation_z_axis[0][2] = 0;                       rotation_z_axis[1][2] = 0;                      rotation_z_axis[2][2] = 1;

//the three vectors defining the x,y,z of "display space" i.e. screen space, used for positioning camera in world space
 	vec d_xvec = mul(rotation_x_axis,mul(rotation_y_axis,mul(rotation_z_axis,vec(1,0,0))));
 	vec d_yvec = mul(rotation_x_axis,mul(rotation_y_axis,mul(rotation_z_axis,vec(0,1,0))));
 	vec d_zvec = mul(rotation_x_axis,mul(rotation_y_axis,mul(rotation_z_axis,vec(0,0,1))));

//this is somewhat abstract - starting from the center of the block, define a sphere
//this sphere has a radius such that the corners of the box are on the surface of the sphere
//multiplying this radius by some amount >1 will move us out from there, a sphere with the same center
//the camera is then located somewhere on this sphere, in a position determined by the use of the rotation matricies

 	//tip radius is the radius of the sphere that touches the tips of the block's corners
 	double tip_radius = std::sqrt(block_xdim_squared+block_ydim_squared+block_zdim_squared)/2.0;
 	double max_dist = 2*2.2*tip_radius;
 	double min_dist = 0.2*tip_radius;
 	//factor of two is for the incremental step of length 0.5
 	//factor of 2.2 is for the tip radius plus the camera sphere radius, 1+1.2

 	vec cam_position = d_center - 1.2*tip_radius*d_yvec; vec cam_up = scale*d_zvec; vec cam_right = scale*d_xvec;

 	int image_center_x = (image_x_dimension-1)/2; int image_current_x;
 	int image_center_y = (image_y_dimension-1)/2; int image_current_y;

 	vec vector_starting_point,vector_test_point;
 	vec vector_increment = 0.5*normalize(-1.0*(cam_position-d_center));

 	vec block_min = vec(0,0,0);
 	vec block_max = vec(block_xdim,block_ydim,block_zdim);

 	Vox temp;
 	RGB	temp_color, curr_color;
  double temp_alpha, curr_alpha;

  std::stack<Vox> empty_voxtack;
  std::stack<Vox> voxtack;

  int alpha_sum;

  unsigned char image_color[3];

 	bool line_box_intersection, color_set;

 	double t0 = 0;
 	double t1 = 9999;

  double tmin, tmax;


 	for(double x = -(image_x_dimension/2-5); x <= (image_x_dimension/2-5); x++)
 		for(double y = -(image_y_dimension/2-5); y <= (image_y_dimension/2-5); y++)
 		{//init (reset)
 			line_box_intersection = false; alpha_sum = 0; color_set = false;    //reset flag values for the new pixel

      curr_alpha = 1.0; curr_color.red = curr_color.green = curr_color.blue = 0; //this is used to process the alpha

      voxtack = empty_voxtack;                                           //reset the stack by setting it equal to an empty stack

 			image_current_x = image_center_x + x; image_current_y = image_center_y + y; //x and y values for the new pixel

 			// if(perspective == true) //this gets added inside the loop - note that the linetest will have to consider the perspective corrected ray
 			// 	vector_increment_perspective = vector_increment + x*0.1*cam_right - y*0.1*cam_up;
      //orthogonal display will have vector increment equal for all pixels

 			vector_starting_point = cam_position + x*cam_up + y*cam_right;

 			//figure out if the parametric line established by parameter z and
 			//	L = vector_starting_point + z*vector_increment
 			//intersects with the box established by (0,0,0) (x,y,z)
 			//i.e. block_min and block_max

 			//The goal here is to know whether or not there is data to sample behind any given pixel - this offers a significant speedup
 			line_box_intersection = parent->intersect_ray_bbox(block_min,block_max,vector_starting_point,vector_increment,tmin,tmax,t0,t1);

 			if(line_box_intersection)
 			{//the ray hits the box, we will need to step through the box
 				for(double z = tmin; z <= tmax; z+=0.5) //go from close intersection point (tmin) to far intersection point (tmax)
 				{
 					vector_test_point = linalg::floor(vector_starting_point + z*vector_increment);  //get the test point
 					temp = parent->get_data_by_vector_index(vector_test_point);                     //get the data at the test point
          alpha_sum += temp.alpha;

          voxtack.push(temp);                                                             //push the data onto the stack

          if(alpha_sum >= 255)
          {
            break;
          }
 				}//end for (z)
        //the for loop is completed, now process the stack

        if(alpha_sum == 0 || voxtack.empty()) //this might be a weird condition to end up in, but just in case
        {
          img.draw_point(image_current_x,image_current_y,black);
          color_set = true;
        }


        if(!color_set)
        {
          while(!voxtack.empty()) //skip this loop if there is no color data, or if the color has already been set
          {//process the stack - math is from https://en.wikipedia.org/wiki/Alpha_compositing
            temp = voxtack.top(); voxtack.pop();
            temp_alpha = temp.alpha / 255;          //map the alpha value to a range between 0 and 1

            curr_color.red = ((parent->palette[temp.state].red * temp_alpha) + (curr_color.red * curr_alpha * (1-temp_alpha)))/(temp_alpha + curr_alpha*(1-temp_alpha));
            curr_color.green = ((parent->palette[temp.state].green * temp_alpha) + (curr_color.green * curr_alpha * (1-temp_alpha)))/(temp_alpha + curr_alpha*(1-temp_alpha));
            curr_color.blue = ((parent->palette[temp.state].blue * temp_alpha) + (curr_color.blue * curr_alpha * (1-temp_alpha)))/(temp_alpha + curr_alpha*(1-temp_alpha));

            curr_alpha = temp_alpha + curr_alpha * (1-temp_alpha);
          }//end while (stack processing)

        //now set the image color based on the computed color
          image_color[0] = curr_color.red;
          image_color[1] = curr_color.green;
          image_color[2] = curr_color.blue;

          img.draw_point(image_current_x,image_current_y,image_color);
        }
 			}
 			else //I saw a ray that did not hit the box, I want to paint it black
 			{
 				img.draw_point(image_current_x,image_current_y,black);
 			}
 		}//end for (x and y)

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

void Voraldo_Draw::init_block(vec dimensions)
{

 if(parent->data != NULL)
 {
		delete[] parent->data;
 }

  parent->x_dim = dimensions[0];
	parent->y_dim = dimensions[1];
	parent->z_dim = dimensions[2];

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

    draw_point(current_point,set,draw,mask);
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

 			draw_line_segment(c1,c2,set,draw,mask);
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
 					draw_point(index,set,draw,mask);
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
 					draw_point(index,set,draw,mask);
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
						draw_point(index,set,draw,mask);
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
            draw_point(index,set,draw,mask);
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
 					draw_point(current,set,draw,mask);
 				}
 			}
 		}
 	}
}



//---------------------------
Voraldo::Voraldo()
{
  io = new Voraldo_IO(this);
  draw = new Voraldo_Draw(this);

  palette = new RGB[256];
  //need to fill in all the data for colors

  //colors

  palette[0] = {0,0,0}; //black - duplicate, but here used to represent 'emtpy'

/*weird desaturated palette "steam lords"
  palette[ 1] = { 33, 59, 37,255};	 //#213b25 dark green
  palette[ 2] = { 58,	96,	74,255}	  //#3a604a medium green
  palette[ 3] = { 79,119, 84,255};	 //#4f7754 light green
  palette[ 4] = {161,159,124,255}; 	//#a19f7c light tan
  palette[ 5] = {119,116,	79,255};	 //#77744f medium tan
  palette[ 6] = {119,	92,	79,255};	 //#775c4f light rose
  palette[ 7] = { 96,	59,	58,255};	 //#603b3a dark rose
  palette[ 8] = { 59,	33,	55,255};	 //#3b2137 purple
  palette[ 9] = { 23,	14,	25,255};	 //#170e19 darkest blue (0)
  palette[10] = { 47,	33,	59,255};	 //#2f213b dark blue (1)
  palette[11] = { 67,	58,	96,255};	 //#433a60 dark blue (2)
  palette[12] = { 79,	82,119,255};	 //#4f5277 dark blue (3)
  palette[13] = {101,115,140,255};	 //#65738c light blue (4)
  palette[14] = {124,148,161,255};	 //#7c94a1 light blue (5)
  palette[15] = {160,185,186,255};	 //#a0b9ba light blue (6)
  palette[16] = {192,209,204,255};	 //#c0d1cc light blue (7)
*/

//REDS

  palette[ 1] = {254,0,0};      //MS Light Red High
  palette[ 2] = {235,138,96};   //MS Light Red Low
  palette[ 3] = {126,0,0};      //MS Dark Red High
  palette[ 4] = {138,54,34};    //MS Dark Red Low
  palette[ 5] = {120,24,38};    //T Dark Red
  palette[ 6] = {165,45,39};    //T Red

//ORANGES

  palette[ 7] = {255,77,0};     //Orange 1
  palette[ 8] = {255,120,30};   //Orange 2
  palette[ 9] = {243,120,43};   //Orange 3
  palette[10] = {201,109,69};   //T Orange

 //YELLOWS

  palette[11] = {255,255,4};    //MS Light Yellow High
  palette[12] = {255,217,63};   //MS Light Yellow Low
  palette[13] = {127,107,0};    //Dark Gold
  palette[14] = {126,126,0};    //MS Dark Yellow High
  palette[15] = {170,92,61};    //MS Dark Yellow Low
  palette[16] = {204,165,98};   //T Yellow Dark
  palette[17] = {207,194,129};  //T Yellow Tan
  palette[18] = {209,202,128};  //T Yellow
  palette[19] = {162,157,107};  //T Tan
  palette[20] = {131,107,63};   //T Light Brown
  palette[21] = {99,73,44};     //T Brown
  palette[22] = {65,51,37};     //T Dark Brown

//GREENS

  palette[23] = {6,255,4};      //MS Light Green High
  palette[24] = {108,217,71};   //MS Light Green Low
  palette[25] = {4,126,0};      //MS Dark Green High
  palette[26] = {12,126,69};    //MS Dark Green Low
  palette[27] = {151,181,138};  //T Light Green
  palette[28] = {101,132,92};   //T Med Green
  palette[29] = {34,58,48};     //T Med-Dark Green
  palette[30] = {32,44,17};     //T Dark Green

//BLUE/INDIGO

  palette[31] = {0,0,126};      //MS Dark Blue High
  palette[32] = {34,52,209};    //MS Dark Blue Low
  palette[33] = {0,0,255};      //MS Light Blue High
  palette[34] = {76,129,251};   //MS Light Blue Low
  palette[35] = {62,62,138};    //T Med Blue
  palette[36] = {76,110,173};   //T Dark blue
  palette[37] = {124,168,213};  //T Blue
  palette[38] = {172,220,241};  //T Light Blue
  palette[39] = {4,126,126};    //MS Dark Teal High
  palette[40] = {68,170,204};   //MS Dark Teal Low
  palette[41] = {6,255,255};    //MS Light Teal High
  palette[42] = {123,226,249};  //MS Light Teal low

//VIOLET

  palette[43] = {254,0,255};    //MS Light Purple High
  palette[44] = {226,61,105};   //MS Light Purple Low
  palette[45] = {82,30,46};     //T Maroon
  palette[46] = {126,0,126};    //MS Dark Purple High
  palette[47] = {92,46,120};    //MS Dark Purple Low
  palette[48] = {88,38,79};     //T Darker Purple
  palette[49] = {80,59,104};    //T Dark Purple
  palette[50] = {133,91,105};   //T Purple
  palette[51] = {223,185,202};  //T Pink

//GREYSCALE

  palette[52] = {255,255,255};  //White
  palette[53] = {190,190,190};  //MS Light Grey High
  palette[54] = {181,181,181};  //MS Light Grey Low
  palette[55] = {126,126,126};  //MS Dark Grey High
  palette[56] = {94,96,110};    //MS Dark Grey Low
  palette[57] = {212,237,237};  //T Lighter Grey
  palette[58] = {134,149,152};  //T Med Light Grey
  palette[59] = {95,99,103};    //T Med Grey
  palette[60] = {58,59,61};     //T Dark-Med Grey
  palette[61] = {40,34,31};     //T Dark Grey
  palette[62] = {0,0,0};        //Black

  data = NULL;//call draw.init_block(x,y,z,noise_fill)
  //to populate the data array
}

Voraldo::~Voraldo()
{
  delete[] io;
  delete[] draw;
  delete[] palette;
  delete[] data;
}

Vox Voraldo::get_data_by_vector_index(vec index)
{
  //std::cout << "beginning" << endl;
  int data_index = index[2]*y_dim*x_dim + index[1]*x_dim + index[0];
  //std::cout << "index calculated" << endl;

  bool x_valid = index[0] < x_dim && index[0] >= 0;
  bool y_valid = index[1] < y_dim && index[1] >= 0;
  bool z_valid = index[2] < z_dim && index[2] >= 0;

  bool point_valid = x_valid && y_valid && z_valid;

  //std::cout << "the index is " << index[0] << " " << index[1] << " " << index[2] << endl;

  //std::cout << "index validated: " << data_index << " versus max array size of " << num_cells << endl;


  if(point_valid)
   return data[data_index];
  else
  {
    Vox default_val;
    default_val.state = 0;
    default_val.alpha = 0;
    default_val.mask = false;
    return default_val;
  }
}

void Voraldo::set_data_by_vector_index(vec index, Vox set, bool draw, bool mask)
{
 int data_index = index[2]*y_dim*x_dim + index[1]*x_dim + index[0];

 bool x_valid = index[0] < x_dim && index[0] >= 0;
 bool y_valid = index[1] < y_dim && index[1] >= 0;
 bool z_valid = index[2] < z_dim && index[2] >= 0;

 bool point_valid = x_valid && y_valid && z_valid;

 if(point_valid)
  if(!data[data_index].mask)
  {
   if(draw)
   {
    data[data_index].state = set.state;
   }
   data[data_index].mask = set.mask;
  }
}

bool Voraldo::planetest(vec plane_point, vec plane_normal, vec test_point)
{
 //return false if the point is above the plane
	//return true if the point is below the plane

	double result = 0.0;

	//equation of plane

	// a (x-x1) + b (y-y1) + c (z-z1) = 0

	double a = plane_normal[0];
	double b = plane_normal[1];
	double c = plane_normal[2];

	double x1 = plane_point[0];
	double y1 = plane_point[1];
	double z1 = plane_point[2];

	double x = test_point[0];
	double y = test_point[1];
	double z = test_point[2];

	result = a * (x - x1) + b * (y - y1) + c * (z - z1);

	return (result < 0)?true:false;
}

bool Voraldo::intersect_ray_bbox(vec bbox_min, vec bbox_max, vec ray_org, vec ray_dir, double &tmin, double &tmax, double t0, double t1)
{/*
 * Ray-box intersection using IEEE numerical properties to ensure that the
 * test is both robust and efficient, as described in:
 *
 *      Amy Williams, Steve Barrus, R. Keith Morley, and Peter Shirley
 *      "An Efficient and Robust Ray-Box Intersection Algorithm"
 *      Journal of graphics tools, 10(1):49-54, 2005
 *
 */
//I pulled this code after three attempts at my own implementation didn't work
  vec bbox[2];
	int sign[3];

	vec inv_direction = vec(1/ray_dir[0],1/ray_dir[1],1/ray_dir[2]);

	sign[0] = (inv_direction[0] < 0);
	sign[1] = (inv_direction[1] < 0);
	sign[2] = (inv_direction[2] < 0);

	bbox[0] = bbox_min;
	bbox[1] = bbox_max;


	//already declared (passed in by reference so that they can be used)
  tmin = (bbox[sign[0]][0] - ray_org[0]) * inv_direction[0];
  tmax = (bbox[1-sign[0]][0] - ray_org[0]) * inv_direction[0];

  double tymin = (bbox[sign[1]][1] - ray_org[1]) * inv_direction[1];
  double tymax = (bbox[1-sign[1]][1] - ray_org[1]) * inv_direction[1];

  if ( (tmin > tymax) || (tymin > tmax) )
    return false;
  if (tymin > tmin)
    tmin = tymin;
  if (tymax < tmax)
    tmax = tymax;

  double tzmin = (bbox[sign[2]][2] - ray_org[2]) * inv_direction[2];
  double tzmax = (bbox[1-sign[2]][2] - ray_org[2]) * inv_direction[2];

  if ( (tmin > tzmax) || (tzmin > tmax) )
    return false;
  if (tzmin > tmin)
    tmin = tzmin;
  if (tzmax < tmax)
    tmax = tzmax;
  return ( (tmin < t1) && (tmax > t0) );

}
