#pragma once
#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>

using namespace Eigen;

MatrixXd F_Null(int n);
MatrixXd Gravity(int n);
MatrixXd reaction_sol(int n, MatrixXd clust, MatrixXd V, float h);