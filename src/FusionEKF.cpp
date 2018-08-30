#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  // cout << "Initializing \n";
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  //the measurement matrix for lidar
  // cout << "Initialize H_laser\n";
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  
  //initialize Jacobian matrix
  // cout << "Initialize Hj\n";
  Hj_ <<  0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0;

  //state covariance matrix P
  // cout << "Initialize P_\n";
	MatrixXd P_ = MatrixXd(4, 4);
	P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 1000;

	//initialize transition matrix F_
  // cout << "Initialize F_\n";
	MatrixXd F_ = MatrixXd(4, 4);
  F_ << 0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0;

  //initialize process noise covariance Q
  // cout << "Initialize Q_\n";
  MatrixXd Q_ = MatrixXd(4, 4);
  Q_ << 0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0;

  //initalize state vector
  // cout << "Initialize x_\n";
  VectorXd x_ = VectorXd(4);
  x_ << 0, 0, 0, 0;

  //Initialize EKF
  ekf_.Init(x_, P_, F_, H_laser_, R_laser_, Q_);
  // cout << "DONE initializing \n";
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    // cout << "EKF: " << '\n';
    // Update the timestamp
    previous_timestamp_ = measurement_pack.timestamp_;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = measurement_pack.raw_measurements_(0);
      double phi = measurement_pack.raw_measurements_(1);
      // float rho_dot = measurement_pack.raw_measurements_(2);

      double px = rho * cos(phi);
      double py = rho * sin(phi);

      /**
       * a radar measurement does not contain enough information 
       * to determine the state variable velocities v_x  v_y
       */
      // cout << "Assign radar X_ \n";
      ekf_.x_ << px, py, 0, 0; 
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      float px = measurement_pack.raw_measurements_(0);
      float py = measurement_pack.raw_measurements_(1);
      ekf_.x_ << px, py, 0, 0;
      // cout << "Assign lidar X_ \n";
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  const double noise_ax = 9;
  const double noise_ay = 9;
  const double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  // cout << "dt: " << dt << '\n';
  const double dt_2 = dt * dt;
  const double dt_3 = dt_2 * dt;
  const double dt_4 = dt_3 * dt;

  previous_timestamp_ = measurement_pack.timestamp_;

  // // cout << "Assign F_ \n";
  ekf_.F_ << 1, 0, dt, 0,
            0, 1, 0, dt,
            0, 0, 1, 0,
            0, 0, 0, 1;
  // // cout << "Assign Q_ \n";
  ekf_.Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0,
            0,  dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay,
            dt_3 / 2 * noise_ax, 0, dt_2 * noise_ax, 0,
            0,  dt_3 / 2 * noise_ay, 0, dt_2 * noise_ay;
  // cout << "Calling Prediction()\n";
  ekf_.Predict(); // Radar and Lidar shares the same state transition func

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    // cout << "Calculating Hj \n";
    Hj_ = tools.CalculateJacobian(ekf_.x_); // Calculate Hj_
    cout << "Hj: " << Hj_ << '\n';
    // cout << "Assign Hj \n";
    ekf_.H_ = Hj_;
    // cout << "Assign R_radar \n";
    ekf_.R_ = R_radar_;
    // cout << "Call UpdateEKF() \n";
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    // cout << "Assign H_laster \n";
    ekf_.H_ = H_laser_;
    // cout << "Assign R_laster \n";
    ekf_.R_ = R_laser_;
    // cout << "Call Update()\n";
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
