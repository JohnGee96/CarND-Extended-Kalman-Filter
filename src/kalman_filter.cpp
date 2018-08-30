#include "kalman_filter.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  MatrixXd Ft = F_.transpose();
 	x_ = F_ * x_;
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::update_(const VectorXd &y,
                           const MatrixXd &H, const MatrixXd &R) {
 
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  MatrixXd Ht = H.transpose();
	MatrixXd S = H * P_ * Ht + R;   // Error Matrix
	MatrixXd Si = S.inverse();
	MatrixXd K = P_ * Ht * Si;      // Kalman Gain

	//new estimate
	x_ = x_ + (K * y);
	P_ = (I - K * H) * P_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred; // Error
  update_(y, H_, R_);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  double rho_dot;
  double px = x_(0);
  double py = x_(1);
  double vx = x_(2);
  double vy = x_(3);

  // Convert to polar coordiate
  double rho = sqrt(px * px + py * py);
  double phi = atan2(py, px);
  if (rho < 0.00001) {
    std::cout << "WARNING - DIVIDE BY ZERO: rho is close to zero\n";
    rho = 0.00001;
  }
  rho_dot = (px * vx + py * vy) / rho;
  
  VectorXd z_pred = VectorXd(3);
  z_pred << rho, phi, rho_dot;

  VectorXd y = z - z_pred;      // Error
  y(1) = atan2(sin(y(1)), cos(y(1)));   // Normalize between [-pi, pi]
  update_(y, H_, R_);
}
