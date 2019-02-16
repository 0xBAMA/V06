#include <iostream>
#include <string>
#include <chrono>
#include "includes/voraldo/v.h"

//stream class shit
using std::cout;
using std::endl;
using std::cin;


//Chrono class aliases
using Clock=std::chrono::high_resolution_clock;
using milliseconds=std::chrono::milliseconds;

int init_x = 257;
int init_y = 257;
int init_z = 257;

vec dimensions(init_x,init_y,init_z);

int plus_x_offset = 75;
int minus_x_offset = -75;

int plus_y_offset = 5;
int minus_y_offset = -5;

int plus_z_offset = 10;
int minus_z_offset = -10;

vec center = vec(floor(init_x/2),floor(init_y/2),floor(init_z/2));

vec a = center + vec(minus_x_offset,plus_y_offset,plus_z_offset);
vec b = center + vec(minus_x_offset,minus_y_offset,plus_z_offset);
vec c = center + vec(plus_x_offset,plus_y_offset,plus_z_offset);
vec d = center + vec(plus_x_offset,minus_y_offset,plus_z_offset);
vec e = center + vec(minus_x_offset,plus_y_offset,minus_z_offset);
vec f = center + vec(minus_x_offset,minus_y_offset,minus_z_offset);
vec g = center + vec(plus_x_offset,plus_y_offset,minus_z_offset);
vec h = center + vec(plus_x_offset,minus_y_offset,minus_z_offset);



Voraldo *main_block;

int main()
{

  main_block = new Voraldo();
  main_block->draw->init_block(dimensions);



  main_block->draw->draw_sphere(center,45,get_vox(37,0.02,false),true,false);

  main_block->draw->draw_quadrilateral_hexahedron(a,b,c,d,e,f,g,h,get_vox(14,0.5,false),true,true);

  main_block->draw->draw_sphere(center,15,get_vox(24,0.14,false),true,true);
  main_block->draw->draw_sphere(center,40,get_vox(0,0,false),true,true);

  main_block->draw->draw_noise(true,0.1);
  //main_block->draw->draw_sphere(vec(floor(init_x/2),floor(init_y/2),floor(init_z/2)+150),20,get_vox(1,255,false),true,false);

  main_block->io->display("new_output.png", 3.14, 3.14/3.0, 3.14/3.0, 0.4, false);

  return 0;
}
