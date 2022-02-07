#include "Auxiliary.h"
#include "Forces_ext.h"


MatrixXd CoM(MatrixXd& V, MatrixXd& M) { 
	//calulus of the Center of Mass
	MatrixXd center(1, 3);
	center.setZero();

	for (int i = 0; i < V.rows(); i++) {
		for (int k = 0; k < 3; k++) {
			center(0, k) += V(i, k) * M(i, 0);
		}
	}
	for (int k = 0; k < 3; k++) {
		center(0, k) /= M.sum();
	}

	return center;
}

void compute_A(MatrixXd& V0, MatrixXd& V1, MatrixXd& M, MatrixXd& A_pq, MatrixXd& A_qq) {

	//centres de masse de l'ancienne et la nouvelle forme
    MatrixXd center0(1, 3);
    MatrixXd center1(1, 3);
    center0 = CoM(V0, M);
    center1 = CoM(V1, M);

    //Calcul de A
    A_pq.setZero(3, 3);
    A_qq.setZero(3, 3);

    for (int i = 0; i < V0.rows(); i++) {
        A_qq += M(i, 0) * (V0.row(i) - center0).transpose() * (V0.row(i) - center0);
        A_pq += M(i, 0) * (V1.row(i) - center1).transpose() * (V0.row(i) - center0);
    }
    A_qq = A_qq.inverse();
}



MatrixXd MatrixSqrt(MatrixXd& S) {
	//renvoit la racine carré d'une matrice
    MatrixXd V(S.rows(), S.rows());
    EigenSolver<MatrixXd> es(S);
    V = es.eigenvectors().real();
    MatrixXd Diag(S.rows(), S.rows());
    Diag = V.inverse() * S * V;

    for (int i = 0; i < S.rows(); i++) {
        Diag(i, i) = std::sqrt(Diag(i, i));
    }
    return V * Diag * V.inverse();
}

void Compute_goals(MatrixXd& A_pq, MatrixXd& A_qq, MatrixXd& V1, MatrixXd& center1, MatrixXd& center2, MatrixXd& G, float beta) {
	//calcul de S
	MatrixXd S(3, 3);
	S = A_pq.transpose() * A_pq;
	S = MatrixSqrt(S);

	//calcul de R
	MatrixXd R(3, 3);
	R = A_pq * S.inverse();

	MatrixXd A(3, 3);
	A = A_pq * A_qq;

	MatrixXd Rbeta = beta * A + (1 - beta) * R;
	//calcul des goals
	MatrixXd goal(1, 3);
	for (int i = 0; i < V1.rows(); i++) {
		goal = (Rbeta * (V1.row(i).transpose() - center1.transpose()) + center2.transpose()).transpose();
		G.row(i) = goal;
	}
}

void Integration(MatrixXd& V_t, MatrixXd& X_t, MatrixXd& goal, MatrixXd& M, double h, double alpha, MatrixXd F_ext, MatrixXd& V_th, MatrixXd& X_th, double damping) {
	//intégration semi-implicite de la vitesse et de la position
	V_th.setZero(V_t.rows(), 3);
	X_th.setZero(X_t.rows(), 3);
	for (int i = 0; i < V_t.rows(); i++) {
		V_th.row(i) = ((1 - damping) * V_t.row(i) + (alpha / h) * (goal.row(i) - X_t.row(i)) + h * F_ext.row(i) / M(i, 0));
		X_th.row(i) = X_t.row(i) + h * V_th.row(i);
	}
}

void transf(MatrixXd Clust1, MatrixXd& Clust2) {
	for (int i = 0; i < Clust1.rows(); i++) {
		MatrixXd a = Clust1.row(i);

		MatrixXd Bend(4, 4);
		Bend.row(0) << 1., 0., 0., 0;
		Bend.row(1) << 0., 1., a(2) * 0.3, 0;
		Bend.row(2) << 0., -a(2) * 0.3, 0.7, 0;
		Bend.row(3) << 0., 0., 0., 1.;

		MatrixXd b(4, 1);
		MatrixXd c(4, 1);

		b << a(0), a(1), a(2), 1.;
		c = Bend * b;
		Clust2.row(i) << c(0), c(1), c(2);
	}
}

