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

int init_x = 256;
int init_y = 128;
int init_z = 128;

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








  vec middle_front_point_1 = vec( 45, 63, 58 ); //side one is towards the negative z
  vec middle_front_point_2 = vec( 45, 63, 69 ); //side one is towards the positive z

  vec middle_back_point_1 = vec( 220, 60, 56 );
  vec middle_back_point_2 = vec( 220, 60, 71 );

  vec outer_point_1 = vec( 210, 32, 38 - 28 );
  vec outer_point_2 = vec( 210, 32, 89 + 28 );

  main_block->draw->draw_triangle( middle_front_point_1, middle_back_point_1, outer_point_1, get_vox(19,1.0,1.0,false));
  main_block->draw->draw_triangle( middle_front_point_2, middle_back_point_2, outer_point_2, get_vox(19,1.0,1.0,false));

  middle_front_point_1 += vec( 1, -1, 0);
  middle_front_point_2 += vec( 1, -1, 0);

  middle_back_point_1 += vec( -1, -1, 0);
  middle_back_point_2 += vec( -1, -1, 0);

  outer_point_1 += vec( -2, -2, 2);
  outer_point_2 += vec( -2, -2, -2);


  main_block->draw->draw_triangle( middle_front_point_1, middle_back_point_1, outer_point_1, get_vox(57,1.0,1.0,false));
  main_block->draw->draw_triangle( middle_front_point_2, middle_back_point_2, outer_point_2, get_vox(57,1.0,1.0,false));


  middle_front_point_1 += vec( 1, -1, 0);
  middle_front_point_2 += vec( 1, -1, 0);

  middle_back_point_1 += vec( -1, -1, 0);
  middle_back_point_2 += vec( -1, -1, 0);

  outer_point_1 += vec( -2, -2, 2);
  outer_point_2 += vec( -2, -2, -2);


  main_block->draw->draw_triangle( middle_front_point_1, middle_back_point_1, outer_point_1, get_vox(61,1.0,1.0,false));
  main_block->draw->draw_triangle( middle_front_point_2, middle_back_point_2, outer_point_2, get_vox(61,1.0,1.0,false));

  middle_front_point_1 += vec( 1, -1, 0);
  middle_front_point_2 += vec( 1, -1, 0);

  middle_back_point_1 += vec( -1, -1, 0);
  middle_back_point_2 += vec( -1, -1, 0);

  outer_point_1 += vec( -2, -2, 2);
  outer_point_2 += vec( -2, -2, -2);


  main_block->draw->draw_triangle( middle_front_point_1, middle_back_point_1, outer_point_1, get_vox(57,1.0,1.0,false));
  main_block->draw->draw_triangle( middle_front_point_2, middle_back_point_2, outer_point_2, get_vox(57,1.0,1.0,false));


  middle_front_point_1 += vec( 1, -1, 0);
  middle_front_point_2 += vec( 1, -1, 0);

  middle_back_point_1 += vec( -1, -1, 0);
  middle_back_point_2 += vec( -1, -1, 0);

  outer_point_1 += vec( -2, -2, 2);
  outer_point_2 += vec( -2, -2, -2);


  main_block->draw->draw_triangle( middle_front_point_1, middle_back_point_1, outer_point_1, get_vox(61,1.0,1.0,false));
  main_block->draw->draw_triangle( middle_front_point_2, middle_back_point_2, outer_point_2, get_vox(61,1.0,1.0,false));
  
  middle_front_point_1 += vec( 1, -1, 0);
  middle_front_point_2 += vec( 1, -1, 0);

  middle_back_point_1 += vec( -1, -1, 0);
  middle_back_point_2 += vec( -1, -1, 0);

  outer_point_1 += vec( -2, -2, 2);
  outer_point_2 += vec( -2, -2, -2);


  main_block->draw->draw_triangle( middle_front_point_1, middle_back_point_1, outer_point_1, get_vox(57,1.0,1.0,false));
  main_block->draw->draw_triangle( middle_front_point_2, middle_back_point_2, outer_point_2, get_vox(57,1.0,1.0,false));


  middle_front_point_1 += vec( 1, -1, 0);
  middle_front_point_2 += vec( 1, -1, 0);

  middle_back_point_1 += vec( -1, -1, 0);
  middle_back_point_2 += vec( -1, -1, 0);

  outer_point_1 += vec( -2, -2, 2);
  outer_point_2 += vec( -2, -2, -2);


  main_block->draw->draw_triangle( middle_front_point_1, middle_back_point_1, outer_point_1, get_vox(61,1.0,1.0,false));
  main_block->draw->draw_triangle( middle_front_point_2, middle_back_point_2, outer_point_2, get_vox(61,1.0,1.0,false));




  //
  // main_block->lighting->apply_ambient_occlusion();
  // main_block->lighting->apply_directional_lighting(5.0, 3.14, 0.25*3.14, 3.14/3, 0.15, true);
  // main_block->lighting->scale_lighting_intensity(4.0);

  main_block->io->save("save.png");

  // main_block->draw->draw_sphere(center, 1000, get_vox(38,0.0007,0.1,false));
  //main_block->draw->draw_maze_base();

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
