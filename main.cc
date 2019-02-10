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
vec center(floor(init_x/2),floor(init_y/2),floor(init_z/2));






Voraldo *main_block;

int main()
{
  Vox vox_type1;
  vox_type1.state = 1;
  vox_type1.alpha = 255;
  vox_type1.mask = false;

  main_block = new Voraldo();
  main_block->draw->init_block(dimensions);
  main_block->draw->draw_sphere(center,5,vox_type1);
  main_block->io->display("new_output.bmp", 3.14, 3.14/3.0, 3.14/3.0, 0.4, false);

  return 0;
}
