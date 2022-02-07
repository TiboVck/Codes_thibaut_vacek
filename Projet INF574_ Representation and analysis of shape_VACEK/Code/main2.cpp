#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>
#include "Auxiliary.h"
#include "Forces_ext.h"

using namespace Eigen; // to use the classes provided by Eigen library

MatrixXd V1; // matrix storing vertex coordinates of the input mesh (n rows, 3 columns) 
MatrixXi F1; // incidence relations between faces and edges (f columns)
MatrixXd V2;
MatrixXi F2;
MatrixXd Clust11;
MatrixXd Clust12;
MatrixXi Fc1;
//MatrixXd Clust21;
//MatrixXd Clust22;
//MatrixXi Fc2;


/**
 * This function is called every time a keyboard button is pressed
 * */
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier) {
  std::cout << "pressed Key: " << key << " " << (unsigned int)key << std::endl;

  if (key == '1') {
    std::cout << "saving to OBJ format" << std::endl;
    igl::writeOBJ("../data/converted_mesh.obj", V1, F1);
  }

  if (key == '2')
  {
    // TO BE COMPLETED

    viewer.data().clear(); // Clear should be called before drawing the mesh
    viewer.data().set_mesh(V1, F1); // update the mesh (both coordinates and faces)
    //viewer.core().align_camera_center(V1, F1);
  }

  return false;
}


void set_meshes(igl::opengl::glfw::Viewer& viewer) {
    viewer.callback_key_down = &key_down; // for dealing with keyboard events
    viewer.data().set_mesh(V1,F1);

    //viewer.append_mesh();
    //viewer.data().set_mesh(Clust22, Fc2);

    
    //viewer.data().set_mesh(V1, F1);

    //viewer.data(0).set_colors(Eigen::RowVector3d(0.3, 0.8, 0.3));
   // viewer.data(1).set_colors(Eigen::RowVector3d(0.8, 0.3, 0.3));

}
/*void set_meshes(igl::opengl::glfw::Viewer& viewer, MatrixXd V, MatrixXi F) {


    viewer.callback_key_down = &key_down; // for dealing with keyboard events
    viewer.data().set_mesh(V, F);

 
    
    //viewer.data().set_mesh(V1, F1);

    //viewer.data(0).set_colors(Eigen::RowVector3d(0.3, 0.8, 0.3));
   // viewer.data(1).set_colors(Eigen::RowVector3d(0.8, 0.3, 0.3));

}*/

// ------------ main program ----------------
int main(int argc, char *argv[])
{//import stars
  igl::readOFF("../data/star.off", V1, F1); // Load an input mesh in OFF format
  igl::readOFF("../data/star_rotated.off", V2, F2);

  //std::cout<<V1<<std::endl;
  //std::cout<<F1<<std::endl;
  
  //V2.row(0) = 1.2 * V2.row(0); //déformation initiale
  //V2.row(1) = 1.6 * V2.row(1);
  //V2.row(2) = 0.9 * V2.row(2);
  int n = 60;
  float hei = 3.;
  int cl = 3;
  int o = 4+(n-1)*4*2;
  float h = 0.1;
  float offset = 0.;

  MatrixXd sol(4,3);
  MatrixXi Fs(2,3);

  sol.row(0)<< -10.,-10.,-0.0;
  sol.row(1)<< -10.,+10.,-0.0;
  sol.row(2)<< +10.,+10.,-0.0;
  sol.row(3)<< +10.,-10.,-0.0;
  Fs.row(0)<<0,1,2;
  Fs.row(1)<<0,2,3;

  
  
  MatrixXd Clust11(n*4,3);
  MatrixXi Fc1(o,3);

  MatrixXd Clust12(n*4,3);

  Clust11.row(0) << 0.0,0.0,offset;
  Clust11.row(1) << 0.0,0.5,offset;
  Clust11.row(2) << 0.5,0.5,offset;
  Clust11.row(3) << 0.5,0.0,offset;
  Fc1.row(0) << 0,1,2;
  Fc1.row(1) << 0,3,2;

  for(int i = 1; i<n;i++){
    Clust11.row(0+i*4) << 0.0,0.0,-i*hei*1./n + offset;
    Clust11.row(1+i*4) << 0.0,0.5,-i*hei*1./n + offset;
    Clust11.row(2+i*4) << 0.5,0.5,-i*hei*1./n + offset;
    Clust11.row(3+i*4) << 0.5,0.0,-i*hei*1./n + offset;

    

    Fc1.row(2+(i-1)*8) << 0+(i-1)*4,1+(i-1)*4,4+(i-1)*4;
    Fc1.row(3+(i-1)*8) << 1+(i-1)*4,5+(i-1)*4,4+(i-1)*4;

    Fc1.row(4+(i-1)*8) << 2+(i-1)*4,1+(i-1)*4,5+(i-1)*4;
    Fc1.row(5+(i-1)*8) << 2+(i-1)*4,5+(i-1)*4,6+(i-1)*4;

    Fc1.row(6+(i-1)*8) << 2+(i-1)*4,3+(i-1)*4,6+(i-1)*4;
    Fc1.row(7+(i-1)*8) << 3+(i-1)*4,6+(i-1)*4,7+(i-1)*4;

    Fc1.row(8+(i-1)*8) << 3+(i-1)*4,0+(i-1)*4,7+(i-1)*4;
    Fc1.row(9+(i-1)*8) << 0+(i-1)*4,4+(i-1)*4,7+(i-1)*4;


}

  int a = n*4-1;
  Fc1.row(o-2) << a,a-1,a-2;
  Fc1.row(o-1) << a,a-2,a-3;


// construction de la matrice de pliage

  transf(Clust11, Clust12);
  //Clust12 = Clust11;


 
  //  print the number of mesh elements
  std::cout << "Vertices: " << V1.rows() << std::endl;
  std::cout << "Faces:    " << F1.rows() << std::endl;

//CLUSTER 1
  MatrixXd Mass(Clust11.rows(), 1);
  Mass.setOnes();
  //Mass(0, 0) = 10;

  //centre de masse
  MatrixXd com1(1, 3);
  com1 = CoM(Clust11, Mass);

  MatrixXd A_pq(3, 3);
  MatrixXd A_qq(3, 3);


// préparer les matrices dont on a besoin
  MatrixXd S(3, 3);
  MatrixXd R(3, 3);
  MatrixXd Vtemp(Clust11.rows(), 3);
  MatrixXd Speedtemp(Clust11.rows(), 3);
  //MatrixXd a_t(Clust11.rows(), 3);

  MatrixXd Speed(Clust11.rows(), 3);
  Speed.setZero();


  MatrixXd G(Clust11.rows(), 3);

  MatrixXd center1(1, 3);
  MatrixXd center2(1, 3);
  center1 = CoM(Clust11, Mass);
  center2 = CoM(Clust12, Mass);
  

  float beta = 0.2;
  int y = n*4/cl;
   
  /*MatrixXd overlap_a(y/4,3);
  MatrixXd speedover_a(y/4,3);
  MatrixXd overlap_b(y/4,3);
  MatrixXd speedover_b(y/4,3);*/
  

  igl::opengl::glfw::Viewer viewer;


  
  int ecart = 0;
  viewer.core().is_animating = true;
  viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer&)->bool // run animation
  {
      

      for(int i=0; i<cl; i++){

        //if(i<cl-1){

          
          
          MatrixXd clust1(y,3);
          MatrixXd clust2(y,3);

            
          
          

          MatrixXd Mass(clust1.rows(), 1);
          Mass.setOnes();

          MatrixXd Vtempi(clust1.rows(), 3);
          MatrixXd Speedtempi(clust1.rows(), 3);
          MatrixXd Speedi(clust1.rows(), 3);
          //MatrixXd a_ti(a_t.rows(), 3);

          MatrixXd Gi(clust1.rows(), 3);
          


          for(int j=i*y;j<i*y+clust1.rows();j++){
            clust1.row(j-i*y) = Clust11.row(j);
            clust2.row(j-i*y) = Clust12.row(j);

            Vtempi.row(j-i*y) = Vtemp.row(j);
            Speedtempi.row(j-i*y) = Speedtemp.row(j);
            Speedi.row(j-i*y) = Speed.row(j);

            /*if(i!=0){
            if((j-i*y)<(y/4)){

              std::cout<<y/4<<"  "<<j-i*y<<std::endl;
              clust2.row(j-i*y) = overlap_a.row(j-i*y);
              Speedi.row(j-i*y) = speedover_a.row(j-i*y);}}

            if(j-i*y>=y){

              std::cout<<y/4<<"  "<<j-i*y<<std::endl;

              overlap_a.row(j-i*y-y) = Clust12.row(j);
              
              speedover_a.row(j-i*y-y) = Speed.row(j);


            }*/
            

          }

          
          

          center1 = CoM(clust1, Mass);
          center2 = CoM(clust2, Mass);

          

          compute_A(clust1, clust2, Mass, A_pq, A_qq);
          

          //compute goals
          Compute_goals(beta , A_pq , A_qq , clust1 ,center1 , center2, Gi);

          MatrixXd F1(clust1.rows() ,3);
          F1.setZero();
          if(ecart!=0){
          

          F1 = reaction_sol(clust1.rows() , clust2 , Speedi,h);
        }
          

          MatrixXd F2 = Gravity(clust1.rows());

    
          
          MatrixXd F = F1+F2;
          F.setZero();
      
          //std::cout<<clust2<<std::endl;
          //intégration
          Integration(Speedi, clust2, Gi, Clust11, Mass, h, 0.1, F, Speedtempi, Vtempi,y*i);

          
          //a_ti = (Speedi-Speedtempi)/h;

          
          //update
          for (int l = 0; l < clust1.rows(); l++) {
              for (int j = 0; j < 3; j++)
              {
                  clust2(l, j) = Vtempi(l, j);
                  Speedi(l, j) = Speedtempi(l, j);
              }}



          
         //if(i==0){
            for(int j=i*y;j<i*y+clust1.rows();j++){
                       Clust12.row(j)=clust2.row(j-i*y);

                      Speed.row(j) = Speedi.row(j-i*y);}

            
        /*  }
          else{
        for(int j=i*y;j<i*y+clust1.rows();j++){


            if(j-i*y<y/4){


              Clust12.row(j) = (overlap_b.row(j-i*y)+clust2.row(j-i*y))/2;


              Speed.row(j) = (speedover_b.row(j-i*y)+clust2.row(j-i*y))/2;}

            else{
              Clust12.row(j)=clust2.row(j-i*y);

              Speed.row(j) = Speedi.row(j-i*y);


            }


          }}*/

       /* for(int u = clust1.rows()-y/4; u < clust1.rows() ;u++){

          std::cout<<y/4<<"  "<<u-(clust1.rows()-y/4)<<std::endl;
              
              overlap_b.row(u-(clust1.rows()-y/4)) = clust2.row(u);
              
              speedover_b.row(u-(clust1.rows()-y/4)) = Speedi.row(u);

          } 



      }

      /*else{

          MatrixXd clust1(y,3);
          MatrixXd clust2(y,3);

            
          
          

          MatrixXd Mass(clust1.rows(), 1);
          Mass.setOnes();

          MatrixXd Vtempi(clust1.rows(), 3);
          MatrixXd Speedtempi(clust1.rows(), 3);
          MatrixXd Speedi(clust1.rows(), 3);
          //MatrixXd a_ti(a_t.rows(), 3);

          MatrixXd Gi(clust1.rows(), 3);
          


          for(int j=i*y;j<i*y+clust1.rows();j++){
            clust1.row(j-i*y) = Clust11.row(j);
            clust2.row(j-i*y) = Clust12.row(j);

          


            Vtempi.row(j-i*y) = Vtemp.row(j);
            Speedtempi.row(j-i*y) = Speedtemp.row(j);
            Speedi.row(j-i*y) = Speed.row(j);
            
            if(j-i*y<y/4){


              clust2.row(j-i*y) = overlap_a.row(j-i*y);
              Speedi.row(j-i*y) = speedover_a.row(j-i*y);}
            

          }

          
          

          center1 = CoM(clust1, Mass);
          center2 = CoM(clust2, Mass);

          

          compute_A(clust1, clust2, Mass, A_pq, A_qq);
          

          //compute goals
          Compute_goals(beta , A_pq , A_qq , clust1 ,center1 , center2, Gi);

          MatrixXd F1(clust1.rows() ,3);
          F1.setZero();
          if(ecart!=0){
          

          F1 = reaction_sol(clust1.rows() , clust2 , Speedi,h);
        }
          

          MatrixXd F2 = Gravity(clust1.rows());

    
          
          MatrixXd F = F1+F2;
          F.setZero();
      
          //std::cout<<clust2<<std::endl;
          //intégration
          Integration(Speedi, clust2, clust1, G, Mass, h, 0.1, F, Speedtempi, Vtempi,y*i);

          
          //a_ti = (Speedi-Speedtempi)/h;

          
          //update
          for (int l = 0; l < clust1.rows(); l++) {
              for (int j = 0; j < 3; j++)
              {
                  clust2(l, j) = Vtempi(l, j);
                  Speedi(l, j) = Speedtempi(l, j);
              }}


          for(int j=i*y;j<i*y+clust1.rows();j++){


            if(j-i*y<y/4){
              Clust12.row(j) = (overlap_b.row(j-i*y)+clust2.row(j-i*y))/2;


              Speed.row(j) = (speedover_b.row(j-i*y)+clust2.row(j-i*y))/2;}

            else{
              Clust12.row(j)=clust2.row(j-i*y);

              Speed.row(j) = Speedi.row(j-i*y);


            }


          }



      }*/




          }

        ecart = 1;


      viewer.data().clear();

      viewer.data(0).set_mesh(Clust12, Fc1);
      viewer.append_mesh();
      viewer.data(1).set_mesh(sol, Fs);
      viewer.data(0).set_colors(Eigen::RowVector3d(0.3, 0.8, 0.3));
      viewer.data(1).set_colors(Eigen::RowVector3d(0.8, 0.3, 0.3));
      //viewer.append_mesh();
      //viewer.data(2).set_mesh(G, Fc1);

    //

      /*viewer.data(0).clear(); // Clear should be called before drawing the mesh
      viewer.data(0).set_mesh(V3, F3);
      viewer.append_mesh();
      viewer.data(0).set_mesh(V1, F1);
      viewer.append_mesh();
      viewer.data(0).set_mesh(V2, F2);*/
      //viewer.data(0).set_colors(Eigen::RowVector3d(0.3, 0.8, 0.3));// update the mesh (both coordinates and faces)}
      return false; 
    };
      
  viewer.launch(); // run the editor
}


/*
  MatrixXd Clust21(8,3);
  Clust21.row(0) = V1.row(0);
  Clust21.row(1) = V1.row(1);
  Clust21.row(2) = V1.row(2);
  Clust21.row(3) = V1.row(3);
  Clust21.row(4) = V1.row(4);
  Clust21.row(5) = V1.row(5);
  Clust21.row(6) = V1.row(6);
  Clust21.row(7) = V1.row(7);

  MatrixXd Clust22(8,3);
  
  Clust22.row(0) = V2.row(0);
  Clust22.row(1) = V2.row(1);
  Clust22.row(2) = V2.row(2);
  Clust22.row(3) = V2.row(3);
  Clust22.row(4) = V2.row(4);
  Clust22.row(5) = V2.row(5);
  Clust22.row(6) = V2.row(6);
  Clust22.row(7) = V2.row(7);

  MatrixXi Fc2(12,3);
  Fc2.row(0) << 0,1,2;
  Fc2.row(1) << 0,2,3;
  Fc2.row(2) << 0,1,4;
  Fc2.row(4) << 1,4,5;
  Fc2.row(5) << 1,2,5;
  Fc2.row(6) << 2,6,5;
  Fc2.row(7) << 2,3,6;
  Fc2.row(8) << 3,6,7;
  Fc2.row(9) << 0,3,4;
  Fc2.row(10) << 3,4,7;
  Fc2.row(11) << 4,5,6;
  Fc2.row(3) << 4,6,7;

    Clust11.row(0) = V1.row(0);
  Clust11.row(1) = V1.row(1);
  Clust11.row(2) = V1.row(2);
  Clust11.row(3) = V1.row(3);
  Clust11.row(4) = V1.row(13);

  MatrixXd Clust12(12,3);
  
  Clust12.row(0) = V2.row(0);
  Clust12.row(1) = V2.row(1);
  Clust12.row(2) = V2.row(2);
  Clust12.row(3) = V2.row(3);
  Clust12.row(4) = V2.row(13);

  MatrixXi Fc1(6,3);
  Fc1.row(0) << 0,1,2;
  Fc1.row(1) << 0,2,3;
  Fc1.row(2) << 0,1,4;
  Fc1.row(3) << 1,2,4;
  Fc1.row(4) << 0,3,4;
  Fc1.row(5) << 3,2,4;*/
