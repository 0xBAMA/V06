#include <iostream>
#include <string>
#include <chrono>
#include <unistd.h>
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

  main_block = new Voraldo();
  main_block->draw->init_block(dimensions);

  vec center = vec(floor(init_x/2),floor(init_y/2),floor(init_z/2));

  auto tick = Clock::now(); //variable to hold start of the timekeeping
  auto tock = Clock::now(); //variable to hold end of timekeeping

   vec ap = vec(  30,  85, 230);
   vec bp = vec(  30,  75, 230);
   vec cp = vec(  75,  85, 230);
   vec dp = vec(  75,  75, 230);
   vec ep = vec(  30,  85,  30);
   vec fp = vec(  30,  75,  30);
   vec gp = vec(  75,  85,  30);
   vec hp = vec(  75,  75,  30);

   vec ar = vec(  75,  85, 230);
   vec br = vec(  75,  75, 230);
   vec cr = vec( 175,  15, 230);
   vec dr = vec( 165,  15, 230);
   vec er = vec(  75,  85, 175);
   vec fr = vec(  75,  75, 175);
   vec gr = vec( 175,  15, 175);
   vec hr = vec( 165,  15, 175);



   vec as,bs,cs,ds,es,fs,gs,hs;
   vec ib,kb;

   vec v1 = vec(2,2,0);
   vec v2 = vec(0,20,0);

   vec ball_position = ((ap + ep + cp + gp) / 4.0) + vec(0,15,0);


   for(int i = 0; i < 315; i += 1)
   {
     std::cout << "frame number " << i;
     tick = Clock::now();

     main_block->draw->draw_sphere(ball_position,10,get_vox(14,1.0,0.3,false),true,false);

     main_block->draw->draw_noise();

     main_block->draw->mask_all_nonzero();

     main_block->draw->draw_sphere(center,1000,get_vox(0,0.0007,0.1,false),true,false);


     main_block->draw->draw_quadrilateral_hexahedron(ap,bp,cp,dp,ep,fp,gp,hp,get_vox(27,1.0,0.3,false),true,false);
     main_block->draw->draw_quadrilateral_hexahedron(ar,br,cr,dr,er,fr,gr,hr,get_vox(27,1.0,0.3,false),true,false);

     for(int it = 0; it < 10; it++)
     {
       as = vec(  75 + it * (10),  85 + it * (-7), 225);
       bs = vec(  75 + it * (10),  78 + it * (-7), 225);
       cs = vec(  85 + it * (10),  85 + it * (-7), 225);
       ds = vec(  85 + it * (10),  78 + it * (-7), 225);
       es = vec(  75 + it * (10),  85 + it * (-7), 180);
       fs = vec(  75 + it * (10),  78 + it * (-7), 180);
       gs = vec(  85 + it * (10),  85 + it * (-7), 180);
       hs = vec(  85 + it * (10),  78 + it * (-7), 180);

       ib = vec(  80 + it * (10),  78 + it * (-7), 227);
       kb = vec(  80 + it * (10),  78 + it * (-7), 178);



       main_block->draw->draw_quadrilateral_hexahedron(as,bs,cs,ds,es,fs,gs,hs,get_vox(9,0.01,1.0,false),true,false);
       main_block->draw->draw_quadrilateral_hexahedron(as+v1,bs+v1,cs+v1,ds+v1,es+v1,fs+v1,gs+v1,hs+v1,get_vox(0,0.01,1.0,false),true,false);

       main_block->draw->draw_cylinder(ib+v2,ib,1.5,get_vox(4,1.0,0.3,false),true,false);
       main_block->draw->draw_sphere(ib+v2,3,get_vox(46,1.0,0.3,false),true,false);
       main_block->draw->draw_sphere(ib+vec(0,5,0),2.5,get_vox(26,1.0,0.3,false),true,false);



       main_block->draw->draw_cylinder(kb+v2,kb,1.5,get_vox(4,1.0,0.3,false),true,false);
       main_block->draw->draw_sphere(kb+v2,3,get_vox(46,1.0,0.3,false),true,false);
       main_block->draw->draw_sphere(kb+vec(0,5,0),2.5,get_vox(26,1.0,0.3,false),true,false);

     }

     main_block->draw->mask_all_nonzero();

     main_block->draw->draw_sphere(center,1000,get_vox(17,0.0007,0.1,false),true,false);


     main_block->draw->draw_line_segment(vec(180,8,227), vec( 80, 78,227), get_vox(30,1.0,1.0,true));
     main_block->draw->draw_line_segment(vec(180,8,178), vec( 80, 78,178), get_vox(30,1.0,1.0,true));

     //main_block->draw->draw_cylinder(center,center+vec(0,20,20),45,get_vox(37,0.02,1.0,false),true,true);
     main_block->draw->draw_tube(center+vec(0,-10,-10),center+vec(0,10,10),45,55,get_vox(37,0.02,0.3,false),true,true);

     main_block->draw->draw_regular_icosahedron(0.01*i,-0.01*i,0,20,center,get_vox(27,1.0,0.3,false/*vertex material*/),2,get_vox(4,1.0,0.3,false/*edge_material*/),1.5,get_vox(37,0.20,0.3,false/*face_material*/));

    // std::cout << 257*257*257*sizeof(Vox) << endl;


     main_block->lighting->apply_ambient_occlusion();
     main_block->lighting->apply_directional_lighting(5.0, 3.14, 0.25*3.14, 3.14/3, 0.15, true);

     main_block->lighting->scale_lighting_intensity(4.0);


      if(i < 10)
      {
        main_block->io->display("animation/new_output00"+ std::to_string(i) +".png",  3.14, 0.01*i*3.14/3.0, 3.14/3.0, 0.4, false);
      }
      else if(i < 100)
      {
        main_block->io->display("animation/new_output0"+ std::to_string(i) +".png",  3.14, 0.01*i*3.14/3.0, 3.14/3.0, 0.4, false);
      }
      else
      {
        main_block->io->display("animation/new_output"+ std::to_string(i) +".png",  3.14, 0.01*i*3.14/3.0, 3.14/3.0, 0.4, false);
      }
      tock = Clock::now();
      std::cout << " took " << std::chrono::duration_cast<milliseconds>(tock-tick).count() << " milliseconds" << endl;

      main_block->draw->mask_unmask_all();
      main_block->draw->clear_all();
   }//end for

  return 0;
}
