#pragma once
#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>
#include "Auxiliary.h"

using namespace Eigen;

void Compute_goals_plasticity(MatrixXd& A_pq, MatrixXd& A_qq, MatrixXd& V1, MatrixXd& center1, MatrixXd& center2, MatrixXd& G, MatrixXd& State, float h, float beta = 0.0, float c_yield = 0.0, float c_creep = 0.0, float c_max = std::numeric_limits<float>::max());
void compute_A_plasticity(MatrixXd& V0, MatrixXd& V1, MatrixXd& M, MatrixXd& A_pq, MatrixXd& A_qq, MatrixXd& State);