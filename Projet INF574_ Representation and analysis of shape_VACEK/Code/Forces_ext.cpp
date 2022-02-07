#include "Forces_ext.h"

MatrixXd F_Null(int n) {
    MatrixXd F(n, 3);
    F.setZero();
    return F;
}

MatrixXd Gravity(int n) {
    MatrixXd F(n, 3);
    F.setZero();
    for (int i = 0; i < n; i++) {
        F(i, 1) = -10;
    }
    return F;
}

MatrixXd reaction_sol(int n, MatrixXd clust, MatrixXd V, float h) {
    MatrixXd F(n, 3);
    F.setZero();
    for (int i = 0; i < n; i++) {
        if (clust(i, 1) <= 0) {
            F(i, 1) = -1.8 * V(i, 1) / h + 10.1;
        }
    }
    return F;
}
