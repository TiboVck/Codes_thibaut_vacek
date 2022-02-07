#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>
#include "Auxiliary.h"
#include "Forces_ext.h"
#include "Quadratic.h"
#include "Plasticity.h"

using namespace Eigen; // to use the classes provided by Eigen library

/**
 * This function is called every time a keyboard button is pressed
 * */

/*
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier) {
  std::cout << "pressed Key: " << key << " " << (unsigned int)key << std::endl;

  if (key == '1') {
    std::cout << "saving to OBJ format" << std::endl;
    //igl::writeOBJ("../data/converted_mesh.obj", V1, F1);
  }

  if (key == '2')
  {
    // TO BE COMPLETED

    viewer.data().clear(); // Clear should be called before drawing the mesh
    viewer.data().set_mesh(V2, F2); // update the mesh (both coordinates and faces)
    //viewer.core().align_camera_center(V1, F1);
  }

  if (key == 'P') {
      selected_point += 1;
      selected_point %= V1.rows();
  }

  if (key == 'M') {
      selected_point -= 1;
      selected_point %= V1.rows();
  }

  if (key == 'U') {
      V2(selected_point, 0) += 0.1;
  }
  if (key == 'D') {
      V2(selected_point, 0) -= 0.1;
  }
  if (key == 'V') {
      V2(selected_point, 1) += 0.1;
  }
  if (key == 'E') {
      V2(selected_point, 0) -= 0.1;
  }
  return false;
}
*/

void set_meshes(igl::opengl::glfw::Viewer& viewer, MatrixXd& V2, MatrixXi& F2) {
    //viewer.callback_key_down = &key_down; // for dealing with keyboard events
    viewer.data().set_mesh(V2, F2);
    //std::cout << V2 << std::endl;
    //viewer.append_mesh();
    //viewer.data(1).set_mesh(V3, F3);
    //viewer.core().align_camera_center(V2 + V3);
    //viewer.data(0).set_colors(Eigen::RowVector3d(0.3, 0.8, 0.3));
    //viewer.data(1).set_colors(Eigen::RowVector3d(0.8, 0.3, 0.3));

}

void linear(MatrixXd& V1, MatrixXi& F1, MatrixXd& V2, MatrixXi& F2) {
    MatrixXd S(3, 3);
    MatrixXd R(3, 3);
    MatrixXd Vtemp(V1.rows(), 3);
    MatrixXd Speedtemp(V1.rows(), 3);
    MatrixXd Speed(V1.rows(), 3);
    MatrixXd Mass(V1.rows(), 1);
    MatrixXd A_pq(3, 3);
    MatrixXd A_qq(3, 3);
    MatrixXd G(V1.rows(), 3);
    MatrixXd center1(1, 3);
    MatrixXd center2(1, 3);
    Speed.setZero();
    Mass.setOnes();
    center1 = CoM(V1, Mass);

    float beta = 0.4;
    double alpha = 0.3;
    double h = 0.01;
    double damping = 0.05;

    igl::opengl::glfw::Viewer viewer;
    set_meshes(viewer, V2, F2);
    viewer.core().is_animating = true;
    viewer.core().camera_view_angle = 30;
    viewer.core().camera_center = Vector3f(2, -15, 0);

    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer&)->bool // run animation
    {
        center2 = CoM(V2, Mass);
        compute_A(V1, V2, Mass, A_pq, A_qq);

        //compute goals
        Compute_goals(A_pq, A_qq, V1, center1, center2, G);
        MatrixXd Force1 = Gravity(V1.rows());
        MatrixXd Force2 = reaction_sol(V1.rows(), V2, Speed, h);
        MatrixXd Force_tot = Force1 + Force2;

        //intégration
        Integration(Speed, V2, G, Mass, h, alpha, Force_tot, Speedtemp, Vtemp, damping);
        //update
        for (int i = 0; i < V1.rows(); i++) {
            for (int j = 0; j < 3; j++)
            {
                V2(i, j) = Vtemp(i, j);
                Speed(i, j) = Speedtemp(i, j);
            }
        }
        viewer.data().clear();
        set_meshes(viewer, V2, F2);
        return false;};
    viewer.launch();
}

void quadratic(MatrixXd& V1, MatrixXi& F1, MatrixXd& V2, MatrixXi& F2) {
    MatrixXd S(3, 3);
    MatrixXd R(3, 3);
    MatrixXd Vtemp(V1.rows(), 3);
    MatrixXd Speedtemp(V1.rows(), 3);
    MatrixXd Speed(V1.rows(), 3);
    MatrixXd Mass(V1.rows(), 1);
    MatrixXd A_pq(3, 3);
    MatrixXd A_qq(3, 3);
    MatrixXd A_pq_quad(3, 9);
    MatrixXd A_qq_quad(3, 9);
    MatrixXd q_quad(V1.rows(), 9);
    MatrixXd G(V1.rows(), 3);
    MatrixXd center1(1, 3);
    MatrixXd center2(1, 3);
    Speed.setZero();
    Mass.setOnes();
    center1 = CoM(V1, Mass);

    float beta = 0.4;
    double alpha = 0.3;
    double h = 0.01;
    double damping = 0.05;

    igl::opengl::glfw::Viewer viewer;
    set_meshes(viewer, V2, F2);
    viewer.core().is_animating = true;
    viewer.core().camera_view_angle = 40;
    viewer.core().camera_center = Vector3f(0, -10, 0);

    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer&)->bool // run animation
    {
        center2 = CoM(V2, Mass);
        compute_q_quad(V1, center1, q_quad);
        compute_A(V1, V2, Mass, A_pq, A_qq);
        compute_A_quad(V1, V2, Mass, q_quad, A_pq_quad, A_qq_quad);
        
        //compute goals
        compute_goals_quad(A_pq_quad, A_qq_quad, q_quad, A_pq, center1, center2, G, beta);
        MatrixXd Force1 = Gravity(V1.rows());
        MatrixXd Force2 = reaction_sol(V1.rows(), V2, Speed, h);
        MatrixXd Force_tot = Force1 + Force2;

        //intégration
        Integration(Speed, V2, G, Mass, h, alpha, Force_tot, Speedtemp, Vtemp, damping);
        //update
        for (int i = 0; i < V1.rows(); i++) {
            for (int j = 0; j < 3; j++){
                V2(i, j) = Vtemp(i, j);
                Speed(i, j) = Speedtemp(i, j);
            }
        }
        viewer.data().clear();
        set_meshes(viewer, V2, F2);
        return false; };
    viewer.launch();
}

void clusters(MatrixXd& V1, MatrixXi& F1, MatrixXd& V2, MatrixXi& F2) {
    int n = 2;
    float hei = 1.;
    int cl = 1;
    int o = 4 + (n - 1) * 4 * 2;
    float h = 0.1;
    float offset = 4.;

    MatrixXd sol(4, 3);
    MatrixXi Fs(2, 3);

    sol.row(0) << -10., -10., 0.0;
    sol.row(1) << -10., +10., 0.0;
    sol.row(2) << +10., +10., 0.0;
    sol.row(3) << +10., -10., 0.0;
    Fs.row(0) << 0, 1, 2;
    Fs.row(1) << 0, 2, 3;

    MatrixXd Clust11(n * 4, 3);
    MatrixXi Fc1(o, 3);
    MatrixXd Clust12(n * 4, 3);

    Clust11.row(0) << 0.0, 0.0, offset;
    Clust11.row(1) << 0.0, 0.5, offset;
    Clust11.row(2) << 0.5, 0.5, offset;
    Clust11.row(3) << 0.5, 0.0, offset;
    Fc1.row(0) << 0, 1, 2;
    Fc1.row(1) << 0, 3, 2;

    for (int i = 1; i < n; i++) {
        Clust11.row(0 + i * 4) << 0.0, 0.0, i* hei * 1. / n + offset;
        Clust11.row(1 + i * 4) << 0.0, 0.5, i* hei * 1. / n + offset;
        Clust11.row(2 + i * 4) << 0.5, 0.5, i* hei * 1. / n + offset;
        Clust11.row(3 + i * 4) << 0.5, 0.0, i* hei * 1. / n + offset;

        Fc1.row(2 + (i - 1) * 8) << 0 + (i - 1) * 4, 1 + (i - 1) * 4, 4 + (i - 1) * 4;
        Fc1.row(3 + (i - 1) * 8) << 1 + (i - 1) * 4, 5 + (i - 1) * 4, 4 + (i - 1) * 4;

        Fc1.row(4 + (i - 1) * 8) << 2 + (i - 1) * 4, 1 + (i - 1) * 4, 5 + (i - 1) * 4;
        Fc1.row(5 + (i - 1) * 8) << 2 + (i - 1) * 4, 5 + (i - 1) * 4, 6 + (i - 1) * 4;

        Fc1.row(6 + (i - 1) * 8) << 2 + (i - 1) * 4, 3 + (i - 1) * 4, 6 + (i - 1) * 4;
        Fc1.row(7 + (i - 1) * 8) << 3 + (i - 1) * 4, 6 + (i - 1) * 4, 7 + (i - 1) * 4;

        Fc1.row(8 + (i - 1) * 8) << 3 + (i - 1) * 4, 0 + (i - 1) * 4, 7 + (i - 1) * 4;
        Fc1.row(9 + (i - 1) * 8) << 0 + (i - 1) * 4, 4 + (i - 1) * 4, 7 + (i - 1) * 4;
    }

    int a = n * 4 - 1;
    Fc1.row(o - 2) << a, a - 1, a - 2;
    Fc1.row(o - 1) << a, a - 2, a - 3;
    // construction de la matrice de pliage
    transf(Clust11, Clust12);

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
    MatrixXd a_t(Clust11.rows(), 3);

    MatrixXd Speed(Clust11.rows(), 3);
    Speed.setZero();

    MatrixXd G(Clust11.rows(), 3);

    MatrixXd center1(1, 3);
    MatrixXd center2(1, 3);
    center1 = CoM(Clust11, Mass);
    center2 = CoM(Clust12, Mass);

    float beta = 0.2;

    igl::opengl::glfw::Viewer viewer;

    int ecart = 0;
    viewer.core().is_animating = true;
    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer&)->bool // run animation
    {
        int y = n * 4 / cl;

        for (int i = 0; i < cl; i++) {
            MatrixXd clust1(y, 3);
            MatrixXd clust2(y, 3);

            MatrixXd Mass(clust1.rows(), 1);
            Mass.setOnes();

            MatrixXd Vtempi(clust1.rows(), 3);
            MatrixXd Speedtempi(clust1.rows(), 3);
            MatrixXd Speedi(clust1.rows(), 3);
            MatrixXd a_ti(a_t.rows(), 3);

            MatrixXd Gi(clust1.rows(), 3);

            for (int j = i * y; j < i * y + clust1.rows(); j++) {
                clust1.row(j - i * y) = Clust11.row(j);
                clust2.row(j - i * y) = Clust12.row(j);

                if (ecart == 0) {
                    clust2(j - i * y, 2) += 0.5 * i;
                    clust1(j - i * y, 2) += 0.5 * i;
                }
                Vtempi.row(j - i * y) = Vtemp.row(j);
                Speedtempi.row(j - i * y) = Speedtemp.row(j);
                Speedi.row(j - i * y) = Speed.row(j);
                a_ti.row(j - i * y) = a_t.row(j);
            }

            center1 = CoM(clust1, Mass);
            center2 = CoM(clust2, Mass);
            compute_A(clust1, clust2, Mass, A_pq, A_qq);
            //compute goals
            Compute_goals(A_pq, A_qq, clust1, center1, center2, Gi, beta);

            MatrixXd F1 = reaction_sol(clust1.rows(), clust2, Speed, h);
            MatrixXd F2 = Gravity(clust1.rows());
            MatrixXd F = F2;
            if (ecart != 0) {
                MatrixXd F = F1;
            }
            F = F_Null(clust1.rows());
            //intégration
            Integration(Speedi, clust2, Gi, Mass, h, 0.2, F, Speedtempi, Vtempi);

            a_ti = (Speedi - Speedtempi) / h;

            //update
            for (int l = 0; l < clust1.rows(); l++) {
                for (int j = 0; j < 3; j++)
                {
                    clust2(l, j) = Vtempi(l, j);
                    Speedi(l, j) = Speedtempi(l, j);
                }
            }

            for (int j = i * y; j < i * y + clust1.rows(); j++) {
                //clust2(j-i*y,2) -= 0.2*i;
                Clust12.row(j) = clust2.row(j - i * y);
                Vtemp.row(j) = Vtempi.row(j - i * y);
                Speedtemp.row(j) = Speedtempi.row(j - i * y);
                Speed.row(j) = Speedi.row(j - i * y);
                a_t.row(j) = a_ti.row(j - i * y);
            }
        }
        ecart = 1;

        viewer.data().clear();

        viewer.data(0).set_mesh(Clust12, Fc1);
        viewer.append_mesh();
        viewer.data(0).set_mesh(Clust11, Fc1);
        //viewer.data(1).set_mesh(sol, Fs);
        return false; };
    viewer.launch();
}

void plasticity(MatrixXd& V1, MatrixXi& F1, MatrixXd& V2, MatrixXi& F2) {
    MatrixXd S(3, 3);
    MatrixXd R(3, 3);
    MatrixXd Vtemp(V1.rows(), 3);
    MatrixXd Speedtemp(V1.rows(), 3);
    MatrixXd Speed(V1.rows(), 3);
    MatrixXd Mass(V1.rows(), 1);
    MatrixXd A_pq(3, 3);
    MatrixXd A_qq(3, 3);
    MatrixXd G(V1.rows(), 3);
    MatrixXd center1(1, 3);
    MatrixXd center2(1, 3);
    MatrixXd State(3, 3);
    State.setIdentity();
    Speed.setZero();
    Mass.setOnes();
    center1 = CoM(V1, Mass);

    float beta = 0;
    double alpha = 0.3;
    double h = 0.001;
    double damping = 0.05;
    float c_yield = 0.5;
    float c_creep = 0.9;
    float c_max = 3;

    igl::opengl::glfw::Viewer viewer;
    set_meshes(viewer, V2, F2);
    viewer.core().is_animating = true;

    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer&)->bool // run animation
    {
        center2 = CoM(V2, Mass);
        compute_A_plasticity(V1, V2, Mass, A_pq, A_qq, State);
        
        //compute goals
        Compute_goals_plasticity(A_pq, A_qq, V1, center1, center2, G, State, h, beta, c_yield, c_creep, c_max);

        MatrixXd F = F_Null(V1.rows());

        //intégration
        Integration(Speed, V2, G, Mass, h, alpha, F, Speedtemp, Vtemp, damping);
        //update
        for (int i = 0; i < V1.rows(); i++) {
            for (int j = 0; j < 3; j++)
            {
                V2(i, j) = Vtemp(i, j);
                Speed(i, j) = Speedtemp(i, j);
            }
        }
        std::cout << State << std::endl;
        viewer.data().clear();
        set_meshes(viewer, V2, F2);
        return false; };
    viewer.launch();
}
// ------------ main program ----------------
int main(int argc, char *argv[])
{
    MatrixXd V1; // matrix storing vertex coordinates of the input mesh (n rows, 3 columns) 
    MatrixXi F1; // incidence relations between faces and edges (f columns)
    MatrixXd V2;
    MatrixXi F2;
    MatrixXd V3;
    MatrixXi F3;
    MatrixXd Clust11;
    MatrixXd Clust12;
    MatrixXi Fc1;

    igl::readOFF("../../data/bunny.off", V1, F1); // Load an input mesh in OFF format
    igl::readOFF("../../data/bunny.off", V2, F2);
    igl::readOFF("../../data/bunny.off", V3, F3);
    for (int i = 0; i < V1.rows(); i++) {
        V2(i, 1) += 0.5;
        V3(i, 1) += 0.5;
        V3(i, 0) += 2;
    }
    quadratic(V1, F1, V2, F2);
    linear(V1, F1, V3, F3);
}
