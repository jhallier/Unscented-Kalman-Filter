#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

const double Pi = 3.141592654;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_x_;

  // initial state vector
  x_ = VectorXd::Zero(n_x_);

  // initial covariance matrix
  P_ = MatrixXd::Identity(n_x_, n_x_);
  P_(0,0) = 0.15;
  P_(1,1) = 0.15;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3; //(10.0/180.0)*Pi;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**  TODO:
  Complete the initialization. See ukf.h for other member properties.
  Hint: one or more values initialized above might be wildly off...
  */

  time_us_ = 0.0;

  // Constant Turn Rate and velocity model, state dimension is 5, augmented state size is 7
  // (x, y, velocity v, yaw rate, yaw dot (first derivative)

  weights_ = VectorXd::Zero(2 * n_aug_ + 1);
  double wi_n = 1.0 / (2*(lambda_+n_aug_));
  weights_.fill(wi_n);
  weights_(0) = lambda_ / (lambda_+n_aug_);
  Xsig_pred_ = MatrixXd::Zero(n_x_, 2*n_aug_ + 1);


}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  // Read in first measurement and set up initial data
  if (!is_initialized_){

	  time_us_ = meas_package.timestamp_;

	  double px, py, v, yaw, yaw_dot;

	  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		  double rho = meas_package.raw_measurements_(0);
		  double phi = meas_package.raw_measurements_(1);
		  double rho_dot = meas_package.raw_measurements_(2);

		  px = cos(phi) * rho;
		  py = sin(phi) * rho;
		  v = 1.0; //initialize speed to about 7km/h - 2m/s
		  yaw = 1.0; //(10.0/180.0)*Pi;
		  yaw_dot = 0.1; // (10.0/180.0)*Pi;
	  }

	  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
		  px = meas_package.raw_measurements_(0);
		  py = meas_package.raw_measurements_(1);
		  v = 1.0;
		  yaw = 1.0;
		  yaw_dot = 0.1;
	  }

	  x_ << px, py, v, yaw, yaw_dot;

	  is_initialized_ = true;
  }

  if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && (use_laser_ == true)){
	  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
	  this->Prediction(delta_t);
	  this->UpdateLidar(meas_package);
	  time_us_ = meas_package.timestamp_;
  }

  if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) && (use_radar_ == true)){
	  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
	  this->Prediction(delta_t);
	  this->UpdateRadar(meas_package);
	  time_us_ = meas_package.timestamp_;
  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

	// Calculate augmented Sigma points
	MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	P_aug(n_aug_-2, n_aug_-2) = std_a_ * std_a_;
	P_aug(n_aug_-1, n_aug_-1) = std_yawdd_ * std_yawdd_;

	MatrixXd P_sqrt = MatrixXd::Zero(n_aug_, n_aug_);
	P_sqrt = P_aug.llt().matrixL(); //Square root of P_aug
	VectorXd x_aug = VectorXd::Zero(n_aug_);
	x_aug.head(n_x_) = x_;

	MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, 2*n_aug_+1);
	Xsig_aug.col(0) = x_aug;
	for (int i=0; i < n_aug_; i++){
		Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * P_sqrt.col(i);
		Xsig_aug.col(i+n_aug_+1) = x_aug - sqrt(lambda_ + n_aug_) * P_sqrt.col(i);
	}

	double dt2 = (delta_t * delta_t) / 2.0;

	//predict sigma points
	for (int i=0; i < (2*n_aug_+1); i++){
		VectorXd x_pred = VectorXd::Zero(n_x_);
		double px = Xsig_aug(0,i);
		double py = Xsig_aug(1,i);
		double v = Xsig_aug(2,i);
		double yaw = Xsig_aug(3,i);
		double yaw_dot = Xsig_aug(4,i);
		double nu_a = Xsig_aug(5,i);
		double nu_phi_2dot = Xsig_aug(6,i);

	  // if yaw_dot zero, use calculation from constant velocity model
		if (abs(yaw_dot) < 0.001){
			x_pred(0) = px + (v * cos(yaw) * delta_t) + (dt2 * cos(yaw) * nu_a);
			x_pred(1) = py + (v * sin(yaw) * delta_t) + (dt2 * sin(yaw) * nu_a);
	    }
	    else{
			x_pred(0) = px + ((v / yaw_dot) * (sin(yaw + yaw_dot * delta_t) - sin(yaw))) + (dt2 * cos(yaw) * nu_a);
			x_pred(1) = py + ((v / yaw_dot) * (-cos(yaw + yaw_dot * delta_t) + cos(yaw))) + (dt2 * sin(yaw) * nu_a);
	    }
		x_pred(2) = v + (delta_t * nu_a);
		x_pred(3) = yaw + (yaw_dot * delta_t) + (dt2 * nu_phi_2dot);
		x_pred(4) = yaw_dot + (delta_t * nu_phi_2dot);

		Xsig_pred_(0,i) = x_pred(0);
		Xsig_pred_(1,i) = x_pred(1);
		Xsig_pred_(2,i) = x_pred(2);
		Xsig_pred_(3,i) = x_pred(3);
		Xsig_pred_(4,i) = x_pred(4);
	}

	//predict state mean
	x_.fill(0.0);
	for (int i=0; i< (2*n_aug_ +1); i++){
		x_ = x_ + weights_(i) * Xsig_pred_.col(i);
	}

	//predict state covariance matrix
	P_.fill(0.0);
	for (int i=0; i< (2*n_aug_+1); i++){
	  VectorXd x_minus = Xsig_pred_.col(i) - x_;
	  // Angle normalization of yaw
	  while (x_minus(3) > Pi) x_minus(3) -= 2.*Pi;
	  while (x_minus(3) < -Pi) x_minus(3) += 2.*Pi;
	  P_ = P_ + weights_(i) * x_minus * x_minus.transpose();
	}

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

	int n_z = 2; // Measurement dimension is 2 (px and py)

	//Transform sigma points into measurement space. Upper two rows are px and py.
	MatrixXd Zsig = MatrixXd::Zero(n_z, 2*n_aug_+1);
	Zsig = Xsig_pred_.topRows(n_z);
	VectorXd z_pred = VectorXd::Zero(n_z);

	// Calculate mean predicted measurement
	z_pred.fill(0.0);
	for (int i=0; i< (2 * n_aug_ + 1); i++){
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	//calculate measurement covariance matrix S
	MatrixXd S = MatrixXd::Zero(n_z, n_z);
	MatrixXd R = MatrixXd::Zero(n_z, n_z);
	R(0,0) = std_laspx_ * std_laspx_;
	R(1,1) = std_laspy_ * std_laspy_;
	// Cross correlation between predicted sigma points and state and measurement mean
	MatrixXd T = MatrixXd::Zero(n_x_, n_z);
	S.fill(0.0);
	T.fill(0.0);

	VectorXd z_delta = VectorXd::Zero(n_z);
	VectorXd x_delta = VectorXd::Zero(n_x_);

	for (int i=0; i < (2*n_aug_+1); i++){
		z_delta = Zsig.col(i) - z_pred;
		x_delta = Xsig_pred_.col(i) - x_;
		while (x_delta(3) > Pi) x_delta(3) -= 2.*Pi;
		while (x_delta(3) < -Pi) x_delta(3) += 2.*Pi;
		T = T + weights_(i) * x_delta * z_delta.transpose();
		S = S + weights_(i) * z_delta * z_delta.transpose();
	}
	S = S + R;

	VectorXd z_meas = meas_package.raw_measurements_;

	// Kalman gain K
	MatrixXd Si = S.inverse();
	MatrixXd K = T * Si;

	// Update state mean and state covariance matrix
	z_delta.fill(0.0);
	z_delta = z_meas - z_pred;

	x_ = x_ + (K * z_delta);
	P_ = P_ - (K * S * K.transpose());

	// LIDAR NIS
	double NIS_lidar = z_delta.transpose() * Si * z_delta;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

	int n_z = 3; //Radar measurement dimension
	MatrixXd Zsig = MatrixXd::Zero(n_z, 2*n_aug_+1);
	VectorXd z_pred = VectorXd::Zero(n_z);

	Zsig.fill(0.0);
	//transform sigma points into measurement space
	for (int i=0; i< (2*n_aug_+1); i++){
	  double px = Xsig_pred_(0, i);
	  double py = Xsig_pred_(1, i);
	  double v = Xsig_pred_(2, i);
	  double yaw = Xsig_pred_(3, i);

	  double rho = sqrt(px*px + py*py);
	  if (abs(rho) < 0.001){
		  std::cout << "rho = " << rho << std::endl;
		  rho = 0.001;
	  }
	  if (abs(px) < 0.001){
		  std::cout << "px = " << px << std::endl;
		  px = 0.001;
	  }
	  double phi = atan2(py,px);
	  double rho_dot = ((px*cos(yaw)*v) + (py*sin(yaw)*v)) / rho;

	  Zsig(0,i) = rho;
	  Zsig(1,i) = phi;
	  Zsig(2,i) = rho_dot;
	}

	//calculate mean predicted measurement
	z_pred.fill(0.0);
	for (int i=0; i< (2*n_aug_ + 1); i++){
	  z_pred = z_pred + weights_(i) * Zsig.col(i);
	}
	while(z_pred(1) > Pi) z_pred(1)-=2*Pi;
	while(z_pred(1) < -Pi) z_pred(1)+=2*Pi;

	//calculate measurement covariance matrix S
	MatrixXd S = MatrixXd::Zero(n_z, n_z);
	MatrixXd R = MatrixXd::Zero(n_z, n_z);
	R.fill(0.0);
	R(0,0) = std_radr_ * std_radr_;
	R(1,1) = std_radphi_ * std_radphi_;
	R(2,2) = std_radrd_ * std_radrd_;
	// Cross correlation between predicted sigma points and state and measurement mean
	MatrixXd T = MatrixXd::Zero(n_x_, n_z);
	T.fill(0.0);
	S.fill(0.0);

	VectorXd z_delta = VectorXd::Zero(n_z);
	VectorXd x_delta = VectorXd::Zero (n_x_);
	for (int i=0; i < (2*n_aug_+1); i++){
		z_delta.fill(0.0);
		z_delta = Zsig.col(i) - z_pred;
		x_delta.fill(0.0);
		x_delta = Xsig_pred_.col(i) - x_;
		while (z_delta(1) > Pi) z_delta(1) -= 2.*Pi;
		while (z_delta(1) < -Pi) z_delta(1) += 2.*Pi;
		while (x_delta(3) > Pi) x_delta(3) -= 2.*Pi;
		while (x_delta(3) < -Pi) x_delta(3) += 2.*Pi;
		T = T + weights_(i) * x_delta * z_delta.transpose();
		S = S + weights_(i) * z_delta * z_delta.transpose();
	}
	S = S + R;

	VectorXd z_meas = VectorXd::Zero(n_z);
	z_meas.fill(0.0);
	z_meas << meas_package.raw_measurements_;

	// Kalman gain K
	MatrixXd K = MatrixXd::Zero(n_x_, n_z);
	MatrixXd Si = S.inverse();
	K = T * Si;

	// Update state mean and state covariance matrix
	z_delta.fill(0.0);
	z_delta = z_meas - z_pred;
	while (z_delta(1) > Pi) z_delta(1) -= 2.*Pi;
	while (z_delta(1) < -Pi) z_delta(1) += 2.*Pi;
	x_ = x_ + (K * z_delta);
	P_ = P_ - (K * S * K.transpose());

	// RADAR NIS
	double NIS_radar = z_delta.transpose() * Si * z_delta;

}
