#pragma once
#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>

using namespace Eigen;

MatrixXd CoM(MatrixXd& V, MatrixXd& M);
void compute_A(MatrixXd& V0, MatrixXd& V1, MatrixXd& M, MatrixXd& A_pq, MatrixXd& A_qq);
MatrixXd MatrixSqrt(MatrixXd& S);
void Integration(MatrixXd& V_t, MatrixXd& X_t, MatrixXd& goal, MatrixXd& M, double h, double alpha, MatrixXd F_ext, MatrixXd& V_th, MatrixXd& X_th, double dampling = 0.0);
void Compute_goals(MatrixXd& A_pq , MatrixXd& A_qq , MatrixXd& V1 , MatrixXd& center1 , MatrixXd& center2, MatrixXd& G, float beta = 0.0);
void transf(MatrixXd Clust1, MatrixXd& Clust2);
