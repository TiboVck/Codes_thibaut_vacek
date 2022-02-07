#include "Plasticity.h"

void Compute_goals_plasticity (MatrixXd& A_pq, MatrixXd& A_qq, MatrixXd& V1, MatrixXd& center1, MatrixXd& center2, MatrixXd& G, MatrixXd& State, float h, float beta, float c_yield, float c_creep, float c_max) {

	//calcul de S
	MatrixXd S(3, 3);
	S = A_pq.transpose() * A_pq;
	S = MatrixSqrt(S);

	MatrixXd I(3, 3);
	I.setIdentity();
	if ((I - State).norm() > c_max) {
		State = (I + c_max * (State - I)) / (State - I).norm();
		State /= std::pow(State.determinant(), 1.0 / 3.0);
	}

	else if ((I - S).norm() > c_yield) {
		State = (I + h * c_creep * (S - I)) * State;
		State /= std::pow(State.determinant(), 1.0 / 3.0);
	}

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
void compute_A_plasticity(MatrixXd& V0, MatrixXd& V1, MatrixXd& M, MatrixXd& A_pq, MatrixXd& A_qq, MatrixXd& State) {
	MatrixXd center0(1, 3);
	MatrixXd center1(1, 3);
	center0 = CoM(V0, M);
	center1 = CoM(V1, M);

	//Calcul de A
	A_pq.setZero(3, 3);
	A_qq.setZero(3, 3);
	MatrixXd p(1, 3);
	MatrixXd q(1, 3);
	for (int i = 0; i < V0.rows(); i++) {
		q = (State * (V0.row(i) - center0).transpose()).transpose();
		p = V1.row(i) - center1;
		A_qq += M(i, 0) * (q.transpose()) * q;
		A_pq += M(i, 0) * (p.transpose()) * q;
	}
	A_qq = A_qq.inverse();
}