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

int init_x = 150;
int init_y = 150;
int init_z = 150;

vec dimensions(init_x,init_y,init_z);

int startval = 0;
int endval = 1;

Voraldo *main_block;

int main()
{
  char numstring[5];

  main_block = new Voraldo();
  main_block->draw->init_block(dimensions);

  vec center = vec(floor(init_x/2),floor(init_y/2),floor(init_z/2));

  auto tick = Clock::now(); //variable to hold start of the timekeeping
  auto tock = Clock::now(); //variable to hold end of timekeeping


  //main_block->io->load_model_from_file("model.mo");

  //std::cout << endl << "loading of model done" << endl;

  //main_block->draw->draw_heightmap();

  // main_block->draw->draw_sphere(center, 1000, get_vox(38,0.0007,0.1,false));
  // main_block->draw->draw_noise();
  //
  // main_block->draw->draw_heightmap();
  //
  // main_block->draw->draw_regular_icosahedron(0,0,0,20,vec(400,100,260),get_vox(27,1.0,0.3,false),2,get_vox(62,0.3,0.3,false),1.8,get_vox(8,0.08,0.3,false));
  //
  // main_block->lighting->apply_ambient_occlusion();
  // main_block->lighting->apply_directional_lighting(5.0, 3.14, 0.25*3.14, 3.14/3, 0.15, true);
  // main_block->lighting->scale_lighting_intensity(4.0);

  main_block->draw->draw_sphere(center, 1000, get_vox(38,0.0007,0.1,false));
  main_block->draw->draw_maze_base();

  for(int i = startval; i < endval; i += 1)
  {
    std::cout << "frame number " << i;
    tick = Clock::now();

    sprintf(numstring, "%04d", i);//this is a better solution on string formatting
    main_block->io->display("animation/new_output"+ std::string(numstring) +".png",  3.14, 0.002*i*3.14/3.0, 3.14/3.0, 0.4, false);

    tock = Clock::now();
    std::cout << " took " << std::chrono::duration_cast<milliseconds>(tock-tick).count() << " milliseconds" << endl;

   }//end for

/*
 if(fork()){
    if(fork()){
      for(int i = startval; i < endval; i += 4)
      {
        std::cout << "frame number " << i;
        tick = Clock::now();

        sprintf(numstring, "%04d", i);//this is a better solution on string formatting
        main_block->io->display("animation/new_output"+ std::string(numstring) +".png",  3.14, 0.002*i*3.14/3.0, 3.14/3.0, 0.4, false);

        tock = Clock::now();
        std::cout << " took " << std::chrono::duration_cast<milliseconds>(tock-tick).count() << " milliseconds" << endl;

       }//end for
   }else{
     for(int i = startval+1; i < endval; i += 4)
     {
       std::cout << "frame number " << i;
       tick = Clock::now();

       sprintf(numstring, "%04d", i);//this is a better solution on string formatting
       main_block->io->display("animation/new_output"+ std::string(numstring) +".png",  3.14, 0.002*i*3.14/3.0, 3.14/3.0, 0.4, false);

       tock = Clock::now();
       std::cout << " took " << std::chrono::duration_cast<milliseconds>(tock-tick).count() << " milliseconds" << endl;

      }//end for
   }
 }else{
   if(fork()){
     for(int i = startval+2; i < endval; i += 4)
     {
       std::cout << "frame number " << i;
       tick = Clock::now();

       sprintf(numstring, "%04d", i);//this is a better solution on string formatting
       main_block->io->display("animation/new_output"+ std::string(numstring) +".png",  3.14, 0.002*i*3.14/3.0, 3.14/3.0, 0.4, false);

       tock = Clock::now();
       std::cout << " took " << std::chrono::duration_cast<milliseconds>(tock-tick).count() << " milliseconds" << endl;

      }//end for
  }else{
    for(int i = startval+3; i < endval; i += 4)
    {
      std::cout << "frame number " << i;
      tick = Clock::now();

      sprintf(numstring, "%04d", i);//this is a better solution on string formatting
      main_block->io->display("animation/new_output"+ std::string(numstring) +".png",  3.14, 0.002*i*3.14/3.0, 3.14/3.0, 0.4, false);

      tock = Clock::now();
      std::cout << " took " << std::chrono::duration_cast<milliseconds>(tock-tick).count() << " milliseconds" << endl;

     }//end for
  }
}
*/

  return 0;
}
