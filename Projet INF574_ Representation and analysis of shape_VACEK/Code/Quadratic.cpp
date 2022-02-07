#include "Quadratic.h"
#include "Auxiliary.h"

void compute_q_quad(MatrixXd& V0, MatrixXd& center0, MatrixXd& q_quad) {
    q_quad.setZero(V0.rows(), 9);
    for (int i = 0; i < V0.rows(); i++) {
        for (int k = 0; k < 3; k++) {
            q_quad(i, k) = (V0(i, k) - center0(0, k));
            q_quad(i, 3 + k) = (V0(i, k) - center0(0, k)) * (V0(i, k) - center0(0, k));
            q_quad(i, 6 + k) = (V0(i, k) - center0(0, k)) * (V0(i, (k + 1) % 3) - center0(0, (k + 1) % 3));
        }
    }
}

void compute_A_quad(MatrixXd& V0, MatrixXd& V1, MatrixXd& M, MatrixXd& q_quad, MatrixXd& A_pq_quad, MatrixXd& A_qq_quad) {
    MatrixXd center1(1, 3);
    A_pq_quad.setZero(3, 9);
    A_qq_quad.setZero(9, 9);
    center1 = CoM(V1, M);
    for (int i = 0; i < V0.rows(); i++) {
        A_qq_quad += M(i, 0) * (q_quad.row(i)).transpose() * (q_quad.row(i));
        A_pq_quad += M(i, 0) * (V1.row(i) - center1).transpose() * (q_quad.row(i));
    }
    A_qq_quad = A_qq_quad.inverse();
}

void compute_goals_quad(MatrixXd& A_pq_quad, MatrixXd& A_qq_quad, MatrixXd& q_quad, MatrixXd& A_pq, MatrixXd& center1, MatrixXd& center2, MatrixXd& G, float beta) {
    
    MatrixXd A(3, 9);
    A = A_pq_quad * A_qq_quad;

    //calcul de S
    MatrixXd S(3, 3);
    S = A_pq.transpose() * A_pq;
    S = MatrixSqrt(S);

    //calcul de R
    MatrixXd R(3, 3);
    R = A_pq * S.inverse();

    //calcul de R
    MatrixXd R_long(3, 9);
    R_long.setZero();
    for (int i = 0; i < 3; i++) {
        R_long.col(i) = R.col(i);
    }

    MatrixXd Rbeta = beta * A + (1 - beta) * R_long;
    //calcul des goals
    MatrixXd goal(1, 3);
    for (int i = 0; i < q_quad.rows(); i++) {
        goal = (Rbeta * (q_quad.row(i)).transpose() + center2.transpose()).transpose();
        G.row(i) = goal;
    }
}
