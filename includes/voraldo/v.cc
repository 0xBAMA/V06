#include "../voraldo/v.h"

using std::cout;
using std::endl;

float clamp(float value, float low, float high)
{
  if(value > high)
  {
    return high;
  }
  else if(value < low)
  {
    return low;
  }
  else
  {
    return value;
  }
}

Vox get_vox(unsigned char state, double alpha, float lighting_intensity, bool mask)
{
  Vox temp;

  temp.state = state;
  temp.alpha = alpha;
  temp.lighting_intensity = lighting_intensity;
  temp.mask = mask;

  return temp;
}


//---------------------------
Voraldo_Lighting::Voraldo_Lighting(Voraldo *p)
{
 parent = p;
}

Voraldo_Lighting::~Voraldo_Lighting()
{

}

void Voraldo_Lighting::apply_directional_lighting(float initial_intensity, double x_rot, double y_rot, double z_rot, double scale, bool divergence)
{
  /*

    herein lies the new directional lighting scheme - it operates very similarly to the display function.
    "divergence" here maps to perspective in the display function

  */

  //cout << "start";
  int image_x_dimension = 3000;
 	int image_y_dimension = 3000;

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

  bool line_box_intersection;
  Vox temp1,temp2;

  double t0 = 0;
  double t1 = 9999;

 double tmin, tmax;

 double alpha_sum; //this holds the value of the light ray

  for(double x = -(image_x_dimension/2-5); x <= (image_x_dimension/2-5); x++)
    for(double y = -(image_y_dimension/2-5); y <= (image_y_dimension/2-5); y++)
    {//init (reset)

      line_box_intersection = false;   //reset flag values for the new pixel

      image_current_x = image_center_x + x; image_current_y = image_center_y + y; //x and y values for the new pixel

      if(divergence == true) //this gets added inside the loop - note that the linetest will have to consider the perspective corrected ray
      {//this isn't working very well
        vector_increment = 0.5*normalize(-1.0*(cam_position-d_center));
        vector_increment = vector_increment + x*0.001*cam_right - y*0.001*cam_up;
      }

      //orthogonal display will have vector_increment equal for all pixels i.e. no divergence
      //figure out if the parametric line established by parameter z and
      //	L = vector_starting_point + z*vector_increment
      //intersects with the box established by (0,0,0) (x,y,z)
      //i.e. block_min and block_max
      vector_starting_point = cam_position + x*cam_up + y*cam_right;


      //The goal here is to know whether or not there is data to sample behind any given pixel - this offers a significant speedup
      line_box_intersection = parent->intersect_ray_bbox(block_min,block_max,vector_starting_point,vector_increment,tmin,tmax,t0,t1);

     if(line_box_intersection)
     {//the ray hits the box, we will need to step through the box
        alpha_sum = 0;
        for(double z = tmin; z <= tmax; z+=0.5) //go from close intersection point (tmin) to far intersection point (tmax)
        {
          vector_test_point = linalg::floor(vector_starting_point + z*vector_increment);  //get the test point
          temp1 = parent->get_data_by_vector_index(vector_test_point);                    //get the data at the test point
          //cout << "test";
          if(temp1.state != 0)
          {
            //set lighting value of the cell, based upon the current lighting intensity, write it back to the array
            //the current lighting intensity is scaled by the input argument initial_intensity, and multiplied by alpha_sum

            //temp2 = get_vox(temp1.state,temp1.alpha,clamp(temp1.lighting_intensity+(initial_intensity*(1-alpha_sum)),0,1),temp1.mask);
            temp2 = get_vox(temp1.state,temp1.alpha,temp1.lighting_intensity+0.001*alpha_sum,temp1.mask);
            parent->set_data_by_vector_index(vector_test_point,temp2,false,false,true);

            //certain amount of light gets 'absorbed' by the cell, temp - subtract the alpha value from the alpha_sum
            alpha_sum += 0.2*temp1.alpha;

            //if the lighting intensity is now less than zero, break out of the for loop - there's no more light left
            if(alpha_sum >= 1.0)
            {
              break;
            }
          }
        }//end for (z)
      }//end if
   }//end for (x and y)
}

void Voraldo_Lighting::apply_ambient_occlusion()
{
  vec index;
  Vox temp;
  double sum;
  int tot;

  for(int i = 0; i < parent->x_dim; i++)
  {
   for(int j = 0; j < parent->y_dim; j++)
   {
     for(int k = 0; k < parent->z_dim; k++)
     {
       sum = 0;
       tot = 0;
       index = vec(i,j,k);


       //this is crude, slow - consider other methods?
       for (int x = -3; x <= 3; x++)
       {
         for (int y = -3; y <= 3; y++)
         {
           for (int z = -3; z <= 3; z++)
           {
             temp = parent->get_data_by_vector_index(vec(i+x,j+y,k+z));
             sum += temp.alpha;
             tot++;
           }
         }
       }
       //do some shit to come up with a number for lighting_intensity
       temp = parent->get_data_by_vector_index(index);
       temp = get_vox(temp.state, temp.alpha, clamp(temp.lighting_intensity*(sum/(tot*2)),temp.lighting_intensity*0.5,1)/*lighting intensity*/, temp.mask);
       parent->set_data_by_vector_index(index,temp,false,false,true);
      }
    }
  }
}


void Voraldo_Lighting::scale_lighting_intensity(double scale)
{
  for(int i = 0; i < parent->num_cells; i++)
  {
   parent->data[i].lighting_intensity =  clamp(parent->data[i].lighting_intensity*scale,0,1);
  }
}

//---------------------------
Voraldo_IO::Voraldo_IO(Voraldo *p)
{
 parent = p;
}

Voraldo_IO::~Voraldo_IO()
{

}

void Voraldo_IO::load_model_from_file(std::string filename)
{//This represents a step in what I think is a good direction - take things out of main and
  //allow for the user to load a model from a list of commands.

  nlohmann::json j,j2; //declare json object
  std::fstream f; //declare fstream object

  int num_steps = 0;
  std::string current_step;
  std::string step_type;

  f.open(filename); //open the user's file
  f >> j; //read in all the content of that file

  num_steps = j.value("num_steps",0); //syntax is (key,default value)

  if(num_steps == 0)
  {
    std::cout << endl << "file says zero steps" << endl;
    return;
  }



  //step 0 - the initialization step
  current_step = j.value("step_0","{\"x\":0,\"y\":0,\"z\":0}");
  j2 = nlohmann::json::parse(current_step);

  parent->draw->init_block(vec(j.value("x",0),j.value("y",0),j.value("z",0)));



  //steps 1 through num_steps - the draw calls
  for(int i = 1; i <= num_steps; i++)
  {
    std::cout << endl << i;
    current_step = j.value("step_" + std::to_string(i),"{\"type\":\"empty\"}");
    j2 = nlohmann::json::parse(current_step);

    step_type = j2.value("type","empty");

    if(step_type == "empty")  //can't switch on a string
    {
      std::cout << endl << "step type is empty" << endl;
    }
    else if(step_type == "unmask_all")
    {
      parent->draw->mask_unmask_all();
    }
    else if(step_type == "invert_mask")
    {
      parent->draw->mask_invert_mask();
    }
    else if(step_type == "mask_all_nonzero")
    {
      parent->draw->mask_all_nonzero();
    }
    else if(step_type == "mask_by_state")
    {
      parent->draw->mask_by_state(j2.value("state",0));
    }
    else if(step_type == "noise")
    {
      parent->draw->draw_noise(j2.value("alpha",0.0),j2.value("lighting_intensity",0.0),j2.value("seed",5),j2.value("draw",true),j2.value("mask",false));
    }
    else if(step_type == "point")
    {
      parent->draw->draw_point(vec(j2.value("point_x",0.0),j2.value("point_y",0.0),j2.value("point_z",0.0)),
                                get_vox(j2.value("set_state",0),j2.value("set_alpha",0.0),j2.value("set_lighting_intensity",0.0),j2.value("set_mask",false)),
                                j2.value("draw",true),j2.value("mask",false));
    }
    else if(step_type == "line_segment")
    {
      parent->draw->draw_line_segment(vec(j2.value("point1_x",0.0),j2.value("point1_y",0.0),j2.value("point1_z",0.0)),
                                      vec(j2.value("point2_x",0.0),j2.value("point2_y",0.0),j2.value("point2_z",0.0)),
                                      get_vox(j2.value("set_state",0),j2.value("set_alpha",0.0),j2.value("set_lighting_intensity",0.0),j2.value("set_mask",false)),
                                      j2.value("draw",true),j2.value("mask",false));
    }
    else if(step_type == "triangle")
    {
      parent->draw->draw_triangle(vec(j2.value("point1_x",0.0),j2.value("point1_y",0.0),j2.value("point1_z",0.0)),
                                  vec(j2.value("point2_x",0.0),j2.value("point2_y",0.0),j2.value("point2_z",0.0)),
                                  vec(j2.value("point3_x",0.0),j2.value("point3_y",0.0),j2.value("point3_z",0.0)),
                                  get_vox(j2.value("set_state",0),j2.value("set_alpha",0.0),j2.value("set_lighting_intensity",0.0),j2.value("set_mask",false)),
                                  j2.value("draw",true),j2.value("mask",false));
    }
    else if(step_type == "sphere")
    {
      parent->draw->draw_sphere(vec(j2.value("point_x",0.0),j2.value("point_y",0.0),j2.value("point_z",0.0)), j2.value("radius",3.0),
                                get_vox(j2.value("set_state",0),j2.value("set_alpha",0.0),j2.value("set_lighting_intensity",0.0),j2.value("set_mask",false)),
                                j2.value("draw",true),j2.value("mask",false));
    }
    else if(step_type == "ellipsoid")
    {
      parent->draw->draw_ellipsoid(vec(j2.value("point_x",0.0),j2.value("point_y",0.0),j2.value("point_z",0.0)),
                                    vec(j2.value("radius_x",0.0),j2.value("radius_y",0.0),j2.value("radius_z",0.0)),
                                    get_vox(j2.value("set_state",0),j2.value("set_alpha",0.0),j2.value("set_lighting_intensity",0.0),j2.value("set_mask",false)),
                                    j2.value("draw",true),j2.value("mask",false));
    }
    else if(step_type == "cylinder")
    {
      parent->draw->draw_cylinder(vec(j2.value("bvec_x",0.0),j2.value("bvec_y",0.0),j2.value("bvec_z",0.0)),
                                  vec(j2.value("tvec_x",0.0),j2.value("tvec_y",0.0),j2.value("tvec_z",0.0)),j2.value("radius",0.0),
                                  get_vox(j2.value("set_state",0),j2.value("set_alpha",0.0),j2.value("set_lighting_intensity",0.0),j2.value("set_mask",false)),
                                  j2.value("draw",true),j2.value("mask",false));
    }
    else if(step_type == "tube")
    {
      parent->draw->draw_tube(vec(j2.value("bvec_x",0.0),j2.value("bvec_y",0.0),j2.value("bvec_z",0.0)),
                                  vec(j2.value("tvec_x",0.0),j2.value("tvec_y",0.0),j2.value("tvec_z",0.0)),j2.value("inner_radius",0.0),j2.value("outer_radius",0.0),
                                  get_vox(j2.value("set_state",0),j2.value("set_alpha",0.0),j2.value("set_lighting_intensity",0.0),j2.value("set_mask",false)),
                                  j2.value("draw",true),j2.value("mask",false));
    }
    else if(step_type == "quadrilateral_hexahedron")
    {
      parent->draw->draw_quadrilateral_hexahedron(vec(j2.value("a_x",0.0),j2.value("a_y",0.0),j2.value("a_z",0.0)),
                                                  vec(j2.value("b_x",0.0),j2.value("b_y",0.0),j2.value("b_z",0.0)),
                                                  vec(j2.value("c_x",0.0),j2.value("c_y",0.0),j2.value("c_z",0.0)),
                                                  vec(j2.value("d_x",0.0),j2.value("d_y",0.0),j2.value("d_z",0.0)),
                                                  vec(j2.value("e_x",0.0),j2.value("e_y",0.0),j2.value("e_z",0.0)),
                                                  vec(j2.value("f_x",0.0),j2.value("f_y",0.0),j2.value("f_z",0.0)),
                                                  vec(j2.value("g_x",0.0),j2.value("g_y",0.0),j2.value("g_z",0.0)),
                                                  vec(j2.value("h_x",0.0),j2.value("h_y",0.0),j2.value("h_z",0.0)),
                                                  get_vox(j2.value("set_state",0),j2.value("set_alpha",0.0),j2.value("set_lighting_intensity",0.0),j2.value("set_mask",false)),
                                                  j2.value("draw",true),j2.value("mask",false));
    }
    else if(step_type == "regular_icosahedron")
    {
      parent->draw->draw_regular_icosahedron(j2.value("x_rot",0.0),j2.value("y_rot",0.0),j2.value("z_rot",0.0),j2.value("scale",0.0),
                                              vec(j2.value("center_x",0.0),j2.value("center_y",0.0),j2.value("center_z",0.0)),
                                              get_vox(j2.value("vertex_state",0),j2.value("vertex_alpha",0.0),j2.value("vertex_lighting_intensity",0.0),j2.value("vertex_mask",false)),j2.value("vertex_radius",0.0),
                                              get_vox(j2.value("edge_state",0),j2.value("edge_alpha",0.0),j2.value("edge_lighting_intensity",0.0),j2.value("edge_mask",false)),j2.value("edge_thickness",0.0),
                                              get_vox(j2.value("face_state",0),j2.value("face_alpha",0.0),j2.value("face_lighting_intensity",0.0),j2.value("face_mask",false)),j2.value("draw_faces",0),
                                              j2.value("draw",true),j2.value("mask",false));
    }
    // else if()
    // {
    //
    // }
    // else if()
    // {
    //
    // }
    // else if()
    // {//if y'all need more
    //
    // }
    else
    {
      std::cout << endl << "unknown draw type" << endl;
    }
    std::cout << " done";
  }
  return;
}

void Voraldo_IO::display(std::string filename, double x_rot, double y_rot, double z_rot, double scale, bool perspective)
{
  using namespace cimg_library;

 	int image_x_dimension = 1919;
 	int image_y_dimension = 1079;

 	CImg<double> img(image_x_dimension,image_y_dimension,1,3,0);

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

  unsigned char image_color[3];

 	bool line_box_intersection, color_set;

 	double t0 = 0;
 	double t1 = 9999;

  double tmin, tmax;
  double alpha_sum;


 	for(double x = -(image_x_dimension/2-5); x <= (image_x_dimension/2-5); x++)
 		for(double y = -(image_y_dimension/2-5); y <= (image_y_dimension/2-5); y++)
 		{//init (reset)
 			line_box_intersection = false; color_set = false;    //reset flag values for the new pixel
      curr_alpha = 1.0; curr_color.red = 0; curr_color.green = 0; curr_color.blue = 0;

      voxtack = empty_voxtack;                             //reset the stack by setting it equal to an empty stack

 			image_current_x = image_center_x + x; image_current_y = image_center_y + y; //x and y values for the new pixel

 			if(perspective == true) //this gets added inside the loop - note that the linetest will have to consider the perspective corrected ray
      {//this isn't working very well
        vector_increment = 0.5*normalize(-1.0*(cam_position-d_center));
 				vector_increment = vector_increment + x*0.001*cam_right - y*0.001*cam_up;
      }
      //orthogonal display will have vector_increment equal for all pixels i.e. no divergence

 			vector_starting_point = cam_position + x*cam_up + y*cam_right;

 			//figure out if the parametric line established by parameter z and
 			//	L = vector_starting_point + z*vector_increment
 			//intersects with the box established by (0,0,0) (x,y,z)
 			//i.e. block_min and block_max

 			//The goal here is to know whether or not there is data to sample behind any given pixel - this offers a significant speedup
 			line_box_intersection = parent->intersect_ray_bbox(block_min,block_max,vector_starting_point,vector_increment,tmin,tmax,t0,t1);

 			if(line_box_intersection)
 			{//the ray hits the box, we will need to step through the box
        alpha_sum = 0;
 				for(double z = tmin; z <= tmax; z+=0.5) //go from close intersection point (tmin) to far intersection point (tmax)
 				{
 					vector_test_point = linalg::floor(vector_starting_point + z*vector_increment);  //get the test point
 					temp = parent->get_data_by_vector_index(vector_test_point);                     //get the data at the test point
          alpha_sum += temp.alpha;
          if(temp.state != 0)
          {
            voxtack.push(temp);
          }
          if(alpha_sum > 4.0) break;
 				}//end for (z)
        //the for loop is completed, now process the stack

        img.draw_point(image_current_x,image_current_y,black);
        while(!voxtack.empty())
        {
          temp = voxtack.top(); voxtack.pop();

          temp_color = parent->palette[temp.state];

          image_color[0] = temp_color.red   * (temp.lighting_intensity);
          image_color[1] = temp_color.green * (temp.lighting_intensity);
          image_color[2] = temp_color.blue  * (temp.lighting_intensity);

          if(temp.alpha != 0.0)
          {
              img.draw_point(image_current_x,image_current_y,image_color,temp.alpha);
          }
        }
 			}
 			else //I saw a ray that did not hit the box, I want to paint it black
 			{
        img.draw_point(image_current_x,image_current_y,black);
 			}
 		}//end for (x and y)

 	img.save_png(filename.c_str());
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
  parent->data[i].alpha = 0.0;
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
  // nlohmann::json j;
  // j["type"] = "unmask_all";
  // std::cout << j.dump() << endl;


 for(int i = 0; i < parent->num_cells; i++)
 {
  parent->data[i].mask = false;
 }
}

void Voraldo_Draw::mask_invert_mask()
{
  // nlohmann::json j;
  // j["type"] = "invert_mask";
  // std::cout << j.dump() << endl;

 for(int i = 0; i < parent->num_cells; i++)
 {
  parent->data[i].mask = parent->data[i].mask ? false:true;
 }
}

void Voraldo_Draw::mask_all_nonzero()
{
  // nlohmann::json j;
  // j["type"] = "mask_all_nonzero";
  // std::cout << j.dump() << endl;

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
  // nlohmann::json j;
  // j["type"] = "mask_by_state";
  // j["state"] = s;
  // std::cout << j.dump() << endl;

 for(int i = 0; i < parent->num_cells; i++)
 {
  if(parent->data[i].state == s)
  {
   parent->data[i].mask = true;
  }
 }
}

void Voraldo_Draw::draw_noise(float alpha, float lighting_intensity, int seed, bool draw, bool mask)
{
  // nlohmann::json j;
  // j["type"] = "noise";
  // j["alpha"] = alpha;
  // j["lighting_intensity"] = lighting_intensity;
  // j["seed"] = seed;
  // j["draw"] = draw;
  // j["mask"] = mask;
  // std::cout << j.dump() << endl;

  std::srand(seed);
  for(int i = 0; i < parent->num_cells; i++)
  {
     if(std::rand()%696 == 69 && std::rand()%696 < 30)
     {
       if(!parent->data[i].mask)
       {
          parent->data[i].state = (std::rand()%30)+1;//this is a little different
          parent->data[i].alpha = alpha;
          parent->data[i].mask = mask;
          parent->data[i].lighting_intensity = lighting_intensity;
       }
     }
  }
}

void Voraldo_Draw::draw_point(vec point, Vox set, bool draw, bool mask)
{
  // nlohmann::json j;
  // j["type"] = "point";
  // j["set_state"] = set.state;
  // j["set_alpha"] = set.alpha;
  // j["set_lighting_intensity"] = set.lighting_intensity;
  // j["set_mask"] = set.mask;
  // j["draw"] = draw;
  // j["mask"] = mask;
  // std::cout << j.dump() << endl;
 parent->set_data_by_vector_index(point,set,draw,mask);
}

void Voraldo_Draw::draw_line_segment(vec v1, vec v2, Vox set, bool draw, bool mask)
{
  // nlohmann::json j;
  // // j["type"] = "line_segment";
  // // j["point1_x"] = v1[0];
  // // j["point1_y"] = v1[1];
  // // j["point1_z"] = v1[2];
  // // j["point2_x"] = v2[0];
  // // j["point2_y"] = v2[1];
  // // j["point2_z"] = v2[2];
  // // j["set_state"] = set.state;
  // // j["set_alpha"] = set.alpha;
  // // j["set_lighting_intensity"] = set.lighting_intensity;
  // // j["set_mask"] = set.mask;
  // // j["draw"] = draw;
  // // j["mask"] = mask;
  // std::cout << j.dump() << endl;

  vec starting_point = v1;
	vec current_point = starting_point;
	vec line_vector = (v2-v1);

	int length = linalg::length(line_vector);

	for(int i = 0; i < length; i++)
	{
		current_point[0] = std::floor(starting_point[0] + i*(line_vector[0]/length));
		current_point[1] = std::floor(starting_point[1] + i*(line_vector[1]/length));
		current_point[2] = std::floor(starting_point[2] + i*(line_vector[2]/length));

    draw_point(current_point,set,draw,mask);
	}
}

void Voraldo_Draw::draw_triangle(vec v0, vec v1, vec v2, Vox set, bool draw, bool mask)
{
  // nlohmann::json j;
  // j["type"] = "triangle";
  // j["point1_x"] = v0[0];
  // j["point1_y"] = v0[1];
  // j["point1_z"] = v0[2];
  // j["point2_x"] = v1[0];
  // j["point2_y"] = v1[1];
  // j["point2_z"] = v1[2];
  // j["point2_x"] = v2[0];
  // j["point2_y"] = v2[1];
  // j["point2_z"] = v2[2];
  // j["set_state"] = set.state;
  // j["set_alpha"] = set.alpha;
  // j["set_lighting_intensity"] = set.lighting_intensity;
  // j["set_mask"] = set.mask;
  // j["draw"] = draw;
  // j["mask"] = mask;
  // std::cout << j.dump() << endl;

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

 		for(double i = 0; i < length; i+=0.3)
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
  // nlohmann::json j;
  // j["type"] = "sphere";
  // j["point_x"] = center_point[0];
  // j["point_y"] = center_point[1];
  // j["point_z"] = center_point[2];
  // j["radius"] = radius;
  // j["set_state"] = set.state;
  // j["set_alpha"] = set.alpha;
  // j["set_lighting_intensity"] = set.lighting_intensity;
  // j["set_mask"] = set.mask;
  // j["draw"] = draw;
  // j["mask"] = mask;
  // std::cout << j.dump() << endl;


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
  // nlohmann::json j;
  // j["type"] = "ellipsoid";
  // j["point_x"] = center_point[0];
  // j["point_y"] = center_point[1];
  // j["point_z"] = center_point[2];
  // j["radius_x"] = radii[0];
  // j["radius_y"] = radii[1];
  // j["radius_z"] = radii[2];
  // j["set_state"] = set.state;
  // j["set_alpha"] = set.alpha;
  // j["set_lighting_intensity"] = set.lighting_intensity;
  // j["set_mask"] = set.mask;
  // j["draw"] = draw;
  // j["mask"] = mask;
  // std::cout << j.dump() << endl;



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
  // nlohmann::json j;
  // j["type"] = "cylinder";
  // j["bvec_x"] = bvec[0];
  // j["bvec_y"] = bvec[1];
  // j["bvec_z"] = bvec[2];
  // j["tvec_x"] = tvec[0];
  // j["tvec_y"] = tvec[1];
  // j["tvec_z"] = tvec[2];
  // j["radius"] = radius;
  // j["set_state"] = set.state;
  // j["set_alpha"] = set.alpha;
  // j["set_lighting_intensity"] = set.lighting_intensity;
  // j["set_mask"] = set.mask;
  // j["draw"] = draw;
  // j["mask"] = mask;
  // std::cout << j.dump() << endl;




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
					if(point_to_line_distance <= radius)
          {
            //cout << "point to line distance test passed" << endl;
            index = vec(i,j,k);
						draw_point(index,set,draw,mask);
					}
				}
			}
		}
	}
}

void Voraldo_Draw::draw_tube(vec bvec, vec tvec, double inner_radius, double outer_radius, Vox set, bool draw, bool mask)
{
  // nlohmann::json j;
  // j["type"] = "tube";
  // j["bvec_x"] = bvec[0];
  // j["bvec_y"] = bvec[1];
  // j["bvec_z"] = bvec[2];
  // j["tvec_x"] = tvec[0];
  // j["tvec_y"] = tvec[1];
  // j["tvec_z"] = tvec[2];
  // j["inner radius"] = inner_radius;
  // j["outer radius"] = outer_radius;
  // j["set_state"] = set.state;
  // j["set_alpha"] = set.alpha;
  // j["set_lighting_intensity"] = set.lighting_intensity;
  // j["set_mask"] = set.mask;
  // j["draw"] = draw;
  // j["mask"] = mask;
  // std::cout << j.dump() << endl;
  //



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

 					if(point_to_line_distance <= outer_radius && point_to_line_distance >= inner_radius)
          {
            index = vec(i,j,k);
            draw_point(index,set,draw,mask);
 					}
 				}
 			}
 		}
 	}
}

void Voraldo_Draw::draw_quadrilateral_hexahedron(vec a, vec b, vec c, vec d, vec e, vec f, vec g, vec h, Vox set, bool draw, bool mask)
{
  // nlohmann::json j;
  // j["type"] = "quadrilateral_hexahedron";
  // j["a_x"] = a[0];
  // j["a_y"] = a[1];
  // j["a_z"] = a[2];
  //
  // j["b_x"] = b[0];
  // j["b_y"] = b[1];
  // j["b_z"] = b[2];
  //
  // j["c_x"] = c[0];
  // j["c_y"] = c[1];
  // j["c_z"] = c[2];
  //
  // j["d_x"] = d[0];
  // j["d_y"] = d[1];
  // j["d_z"] = d[2];
  //
  //
  // j["e_x"] = e[0];
  // j["e_y"] = e[1];
  // j["e_z"] = e[2];
  //
  // j["f_x"] = f[0];
  // j["f_y"] = f[1];
  // j["f_z"] = f[2];
  //
  // j["g_x"] = g[0];
  // j["g_y"] = g[1];
  // j["g_z"] = g[2];
  //
  // j["h_x"] = h[0];
  // j["h_y"] = h[1];
  // j["h_z"] = h[2];
  //
  // j["set_state"] = set.state;
  // j["set_alpha"] = set.alpha;
  // j["set_lighting_intensity"] = set.lighting_intensity;
  // j["set_mask"] = set.mask;
  // j["draw"] = draw;
  // j["mask"] = mask;
  // std::cout << j.dump() << endl;





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

void Voraldo_Draw::draw_regular_icosahedron(double x_rot, double y_rot, double z_rot, double scale, vec center_point, Vox vertex_material, double verticies_radius, Vox edge_material, double edge_thickness, Vox face_material,bool draw_faces, bool draw, bool mask)
{
  // nlohmann::json j2;
  // j2["type"] = "regular_icosahedron";
  //
  // j2["x_rot"] = x_rot;
  // j2["y_rot"] = y_rot;
  // j2["z_rot"] = z_rot;
  // j2["scale"] = scale;
  // j2["center_x"] = center_point[0];
  // j2["center_y"] = center_point[1];
  // j2["center_z"] = center_point[2];
  // j2["vertex_state"] = vertex_material.state;
  // j2["vertex_alpha"] = vertex_material.alpha;
  // j2["vertex_lighting_intensity"] = vertex_material.lighting_intensity;
  // j2["vertex_mask"] = vertex_material.mask;
  // j2["vertex_radius"] = verticies_radius;
  // j2["edge_state"] = edge_material.state;
  // j2["edge_alpha"] = edge_material.alpha;
  // j2["edge_lighting_intensity"] = edge_material.lighting_intensity;
  // j2["edge_mask"] = edge_material.mask;
  // j2["edge_thickness"] = edge_thickness;
  // j2["face_state"] = face_material.state;
  // j2["face_alpha"] = face_material.alpha;
  // j2["face_lighting_intensity"] = face_material.lighting_intensity;
  // j2["face_mask"] = face_material.mask;
  // j2["draw_faces"] = draw_faces;
  //
  // j2["draw"] = draw;
  // j2["mask"] = mask;
  // std::cout << j2.dump() << endl;



  double phi = (1 + std::sqrt(5.0))/2.0;

//rotation matricies allowing rotation of the polyhedron
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

  mat rotation = mul(rotation_x_axis,mul(rotation_y_axis,rotation_z_axis)); //multiply these all together

  vec a,b,c,d,e,f,g,h,i,j,k,l;
//the work for this is in my journal, the entry on 3/1/2019
//it is based on the idea that the points of a regular icosahedron lie on the points defined by three mutually orthogonal golden rectangles that share a center point at the origin
//i.e. these rectanges are abcd,efgh and ijkl
  a = mul(rotation, vec(  0,  1*scale,  phi*scale)) + center_point; e = mul(rotation, vec(  1*scale,  phi*scale,  0)) + center_point; i = mul(rotation, vec(  phi*scale,  0,  1*scale)) + center_point;
  b = mul(rotation, vec(  0,  1*scale, -phi*scale)) + center_point; f = mul(rotation, vec( -1*scale, -phi*scale,  0)) + center_point; j = mul(rotation, vec(  phi*scale,  0, -1*scale)) + center_point;
  c = mul(rotation, vec(  0, -1*scale,  phi*scale)) + center_point; g = mul(rotation, vec( -1*scale,  phi*scale,  0)) + center_point; k = mul(rotation, vec( -phi*scale,  0,  1*scale)) + center_point;
  d = mul(rotation, vec(  0, -1*scale, -phi*scale)) + center_point; h = mul(rotation, vec(  1*scale, -phi*scale,  0)) + center_point; l = mul(rotation, vec( -phi*scale,  0, -1*scale)) + center_point;
//nonzero components of the coordinates are scaled by the scale input argument. The result of that operation is multiplied by the composed rotation matrix, then added to the shape's center point


  if(draw_faces)
  {//draw the faces -
    draw_triangle( a, g, e, face_material, draw, mask); //AGE
    draw_triangle( a, i, e, face_material, draw, mask); //AIE
    draw_triangle( a, c, i, face_material, draw, mask); //ACI
    draw_triangle( a, c, k, face_material, draw, mask); //ACK
    draw_triangle( a, g, k, face_material, draw, mask); //AGK
    draw_triangle( l, b, g, face_material, draw, mask); //LBG
    draw_triangle( l, g, k, face_material, draw, mask); //LGK
    draw_triangle( l, f, k, face_material, draw, mask); //LFK
    draw_triangle( l, d, f, face_material, draw, mask); //LDF
    draw_triangle( l, d, b, face_material, draw, mask); //LDB
    draw_triangle( k, f, c, face_material, draw, mask); //KFC
    draw_triangle( f, h, c, face_material, draw, mask); //FHC
    draw_triangle( h, i, c, face_material, draw, mask); //HIC
    draw_triangle( e, j, i, face_material, draw, mask); //EJI
    draw_triangle( b, g, e, face_material, draw, mask); //BGE
    draw_triangle( f, h, d, face_material, draw, mask); //FHD
    draw_triangle( d, h, j, face_material, draw, mask); //DHJ
    draw_triangle( d, b, j, face_material, draw, mask); //DBJ
    draw_triangle( b, j, e, face_material, draw, mask); //BJE
    draw_triangle( h, i, j, face_material, draw, mask); //HIJ
  }




  if(edge_thickness)
  {//nonzero value passed for edge thickness
    //draw the edges
    if(edge_thickness <= 1)
    {//use lines
      draw_line_segment(a,c, edge_material, draw, mask); //AC
      draw_line_segment(a,e, edge_material, draw, mask); //AE
      draw_line_segment(a,g, edge_material, draw, mask); //AG
      draw_line_segment(a,i, edge_material, draw, mask); //AI
      draw_line_segment(a,k, edge_material, draw, mask); //AK

      draw_line_segment(b,d, edge_material, draw, mask); //BD
      draw_line_segment(b,e, edge_material, draw, mask); //BE
      draw_line_segment(b,g, edge_material, draw, mask); //BG
      draw_line_segment(b,j, edge_material, draw, mask); //BJ
      draw_line_segment(b,l, edge_material, draw, mask); //BL

      draw_line_segment(c,f, edge_material, draw, mask); //CF
      draw_line_segment(c,h, edge_material, draw, mask); //CH
      draw_line_segment(c,i, edge_material, draw, mask); //CI
      draw_line_segment(c,k, edge_material, draw, mask); //CK

      draw_line_segment(d,f, edge_material, draw, mask); //DF
      draw_line_segment(d,h, edge_material, draw, mask); //DH
      draw_line_segment(d,j, edge_material, draw, mask); //DJ
      draw_line_segment(d,l, edge_material, draw, mask); //DL

      draw_line_segment(e,g, edge_material, draw, mask); //EG
      draw_line_segment(e,i, edge_material, draw, mask); //EI
      draw_line_segment(e,j, edge_material, draw, mask); //EJ

      draw_line_segment(f,h, edge_material, draw, mask); //FH
      draw_line_segment(f,k, edge_material, draw, mask); //FK
      draw_line_segment(f,l, edge_material, draw, mask); //FL

      draw_line_segment(g,k, edge_material, draw, mask); //GK
      draw_line_segment(g,l, edge_material, draw, mask); //GL

      draw_line_segment(h,i, edge_material, draw, mask); //HI
      draw_line_segment(h,j, edge_material, draw, mask); //HJ

      draw_line_segment(i,j, edge_material, draw, mask); //IJ

      draw_line_segment(k,l, edge_material, draw, mask); //KL
    }
    else
    {//use cylinders of radius edge_thickness
      draw_cylinder(a,c, edge_thickness, edge_material, draw, mask); //AC
      draw_cylinder(a,e, edge_thickness, edge_material, draw, mask); //AE
      draw_cylinder(a,g, edge_thickness, edge_material, draw, mask); //AG
      draw_cylinder(a,i, edge_thickness, edge_material, draw, mask); //AI
      draw_cylinder(a,k, edge_thickness, edge_material, draw, mask); //AK

      draw_cylinder(b,d, edge_thickness, edge_material, draw, mask); //BD
      draw_cylinder(b,e, edge_thickness, edge_material, draw, mask); //BE
      draw_cylinder(b,g, edge_thickness, edge_material, draw, mask); //BG
      draw_cylinder(b,j, edge_thickness, edge_material, draw, mask); //BJ
      draw_cylinder(b,l, edge_thickness, edge_material, draw, mask); //BL

      draw_cylinder(c,f, edge_thickness, edge_material, draw, mask); //CF
      draw_cylinder(c,h, edge_thickness, edge_material, draw, mask); //CH
      draw_cylinder(c,i, edge_thickness, edge_material, draw, mask); //CI
      draw_cylinder(c,k, edge_thickness, edge_material, draw, mask); //CK

      draw_cylinder(d,f, edge_thickness, edge_material, draw, mask); //DF
      draw_cylinder(d,h, edge_thickness, edge_material, draw, mask); //DH
      draw_cylinder(d,j, edge_thickness, edge_material, draw, mask); //DJ
      draw_cylinder(d,l, edge_thickness, edge_material, draw, mask); //DL

      draw_cylinder(e,g, edge_thickness, edge_material, draw, mask); //EG
      draw_cylinder(e,i, edge_thickness, edge_material, draw, mask); //EI
      draw_cylinder(e,j, edge_thickness, edge_material, draw, mask); //EJ

      draw_cylinder(f,h, edge_thickness, edge_material, draw, mask); //FH
      draw_cylinder(f,k, edge_thickness, edge_material, draw, mask); //FK
      draw_cylinder(f,l, edge_thickness, edge_material, draw, mask); //FL

      draw_cylinder(g,k, edge_thickness, edge_material, draw, mask); //GK
      draw_cylinder(g,l, edge_thickness, edge_material, draw, mask); //GL

      draw_cylinder(h,i, edge_thickness, edge_material, draw, mask); //HI
      draw_cylinder(h,j, edge_thickness, edge_material, draw, mask); //HJ

      draw_cylinder(i,j, edge_thickness, edge_material, draw, mask); //IJ

      draw_cylinder(k,l, edge_thickness, edge_material, draw, mask); //KL
    }
  }

  if(verticies_radius)
  {//nonzero value passed for vertex radius
    //draw the vertexes
    draw_sphere(a, verticies_radius, vertex_material, draw, mask);
    draw_sphere(b, verticies_radius, vertex_material, draw, mask);
    draw_sphere(c, verticies_radius, vertex_material, draw, mask);
    draw_sphere(d, verticies_radius, vertex_material, draw, mask);

    draw_sphere(e, verticies_radius, vertex_material, draw, mask);
    draw_sphere(f, verticies_radius, vertex_material, draw, mask);
    draw_sphere(g, verticies_radius, vertex_material, draw, mask);
    draw_sphere(h, verticies_radius, vertex_material, draw, mask);

    draw_sphere(i, verticies_radius, vertex_material, draw, mask);
    draw_sphere(j, verticies_radius, vertex_material, draw, mask);
    draw_sphere(k, verticies_radius, vertex_material, draw, mask);
    draw_sphere(l, verticies_radius, vertex_material, draw, mask);
  }
}

void Voraldo_Draw::draw_heightmap(/*std::string filename, std::vector<Vox> materials,*/ bool draw, bool mask)
{
  using namespace cimg_library;
  unsigned char current_color;
  unsigned char current_height;
  unsigned char current_colormap;

  CImg<unsigned char> heightmap("heights.png");
  CImg<unsigned char> colormap("greycolors.png");

  for(int x = 0; x < 512; x++){
    for(int z = 0; z < 512; z++){

      current_height = heightmap.atXY(x,z);
      current_colormap = colormap.atXY(x,z);

      if(current_colormap < 30){
        current_color = 62; //black
      }else if(current_colormap < 60){
        current_color = 61; // t dark grey
      }else if(current_colormap < 90){
        current_color = 60; // t dark med
      }else if(current_colormap < 120){
        current_color = 59; // t med light
      }else if(current_colormap < 150){
        current_color = 58; // t light
      }else{
        current_color = 57; // t lighter
      }

      draw_point(vec(x,current_height/3,z),get_vox(current_color,1,0.3,false),draw,mask);


      //for(int y = 0; y < current_height/3; y++){
        //draw_point(vec(x,y,z),get_vox(current_color,1,0.3,false),draw,mask);
      //}
      //cout << heightmap.atXY(x,z) << " ";
    }
    //cout << endl;
  }


    //go up to y

  return;
}


void draw_maze_base()
{
  std::vector<CImg<unsigned char>> levelvec;
  //allows the use of a loop

  CImg<unsigned char> level0("maze/level0.png");
  levelvec.push_back(level0);

  CImg<unsigned char> level1("maze/level1.png");
  levelvec.push_back(level1);

  CImg<unsigned char> level2("maze/level2.png");
  levelvec.push_back(level2);

  CImg<unsigned char> level3("maze/level3.png");
  levelvec.push_back(level3);

  CImg<unsigned char> level4("maze/level4.png");
  levelvec.push_back(level4);

  CImg<unsigned char> level5("maze/level5.png");
  levelvec.push_back(level5);

  Vox wall_mat    = get_vox(59,0.3,1.0,false);
  Vox column_mat  = get_vox(13,0.5,1.0,false);
  Vox up_mat      = get_vox(24,1.0,1.0,false);
  Vox down_mat    = get_vox( 5,1.0,1.0,false);

  int z = 25;
  unsigned char red, green, blue;

  for(int n = 0; n < 6; n++)
  {//0,1,2,3,4,5
    for(int x = 0; x < 100; x++)
    {
      for(int y = 0; y < 100; y++)
      {

        red   = levelvec[n].atXY(x,y,0,0);//accessing separate channels
        green = levelvec[n].atXY(x,y,0,1);
        blue  = levelvec[n].atXY(x,y,0,2);


        if(red == 0 && green == 0 && blue == 0)
        {//if it's a black pixel, draw the wall

        }
        else if(red == 255 && green == 255 && blue == 255)
        {//if it's a white pixel, it's negative space, skip it

        }
        else if(red == 255 && green == 0 && blue == 0)
        {//if it's a red pixel, it's a way to go down, draw a column down in the negative space

        }
        else if(red == 0 && green == 255 && blue == 0)
        {//if it's a green pixel, it's a way to go up

        }
      }
    }
    z+=18;
  }
}


//---------------------------
Voraldo::Voraldo()
{
  io = new Voraldo_IO(this);
  draw = new Voraldo_Draw(this);
  lighting = new Voraldo_Lighting(this);

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
    default_val.alpha = 0.0;
    default_val.mask = false;
    return default_val;
  }
}

void Voraldo::set_data_by_vector_index(vec index, Vox set, bool draw, bool mask, bool force)
{
 int data_index = index[2]*y_dim*x_dim + index[1]*x_dim + index[0];

 bool x_valid = index[0] < x_dim && index[0] >= 0;
 bool y_valid = index[1] < y_dim && index[1] >= 0;
 bool z_valid = index[2] < z_dim && index[2] >= 0;

 bool point_valid = x_valid && y_valid && z_valid;

 if(point_valid)
 {
  if(!data[data_index].mask)
  {
   if(draw)
   {
    data[data_index] = set;
   }
    data[data_index].mask = mask; //this takes precedence over the Vox value of mask
  }

  if(force)
  {
    data[data_index] = set;
  }
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
