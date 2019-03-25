#include <iostream>
#include <string>
#include <chrono>
#include <unistd.h>
#include <cstdio>
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

Voraldo *main_block;

int main()
{
  char * numstring;

  main_block = new Voraldo();
  main_block->draw->init_block(dimensions);

  vec center = vec(floor(init_x/2),floor(init_y/2),floor(init_z/2));

  auto tick = Clock::now(); //variable to hold start of the timekeeping
  auto tock = Clock::now(); //variable to hold end of timekeeping


  main_block->io->load_model_from_file("model.mo");

  main_block->lighting->apply_ambient_occlusion();
  main_block->lighting->apply_directional_lighting(5.0, 3.14, 0.25*3.14, 3.14/3, 0.15, true);
  main_block->lighting->scale_lighting_intensity(4.0);

  for(int i = 0; i < 1; i += 1)
  {
    std::cout << "frame number " << i;
    tick = Clock::now();

    sprintf(numstring, "%04d", i);//this is a better solution on string formatting
    main_block->io->display("animation/new_output"+ std::string(numstring) +".png",  3.14, 0.002*i*3.14/3.0, 3.14/3.0+0.02*std::sin(0.1*i), 0.4, false);

    tock = Clock::now();
    std::cout << " took " << std::chrono::duration_cast<milliseconds>(tock-tick).count() << " milliseconds" << endl;

   }//end for

  return 0;
}
