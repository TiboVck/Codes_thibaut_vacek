#pragma once
#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>

using namespace Eigen;

void compute_q_quad(MatrixXd& V0, MatrixXd& center0, MatrixXd& q_quad);
void compute_A_quad(MatrixXd& V0, MatrixXd& V1, MatrixXd& M, MatrixXd& q_quad, MatrixXd& A_pq_quad, MatrixXd& A_qq_quad);
void compute_goals_quad(MatrixXd& A_pq_quad, MatrixXd& A_qq_quad, MatrixXd& q_quad, MatrixXd& A_pq, MatrixXd& center1, MatrixXd& center2, MatrixXd& G, float beta = 0);
