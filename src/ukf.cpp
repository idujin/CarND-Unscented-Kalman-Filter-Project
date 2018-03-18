#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */

UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.7;//30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;//30;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  weights_ = ComputeWeights();

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
  if (!is_initialized_) {

    x_ << 0, 0, 0, 0, 0;
    P_ << 
      1, 0, 0, 0, 0,
      0, 1, 0, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 0, 0, 1;

    float px = 0.f;
    float py = 0.f;
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      float rho = meas_package.raw_measurements_(0);
      float phi = meas_package.raw_measurements_(1);
      px = rho * cosf(phi);
      py = rho * sinf(phi);
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      px = meas_package.raw_measurements_(0);
      py = meas_package.raw_measurements_(1);
    }
    
    x_(0) = px;
    x_(1) = py;
    time_us_ = meas_package.timestamp_;

    is_initialized_ = true;
    return;

  }

  double dt = (meas_package.timestamp_ - time_us_) / 1000000.f;
  time_us_ = meas_package.timestamp_ ;

  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_  + 1);

  Prediction(&Xsig_pred, &x_, &P_, weights_, std_a_, std_yawdd_, dt);
  PredictMeanAndCovariance(&x_, &P_, Xsig_pred, weights_);

  if(meas_package.sensor_type_ == MeasurementPackage::RADAR){

    UpdateRadar(&x_,
                &P_,
                meas_package,
                Xsig_pred,
                weights_,
                std_radr_,
                std_radphi_,
                std_radrd_);

  }
  else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
    UpdateLidar(&x_,
                &P_,
                meas_package,
                Xsig_pred,
                weights_,
                std_laspx_,
                std_laspy_);
  }

  //cout << "x_ = " << x_ << endl;
  //cout << "P_ = " << P_ << endl;
  
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction( MatrixXd* Xsig_pred_out,
                      VectorXd* x,
                      MatrixXd* P,
                      VectorXd weights,
                      const double std_a, 
                      const double std_yawdd,
                      double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  VectorXd x_in = *x;
  // initial covariance matrix
  MatrixXd P_in = *P;

  MatrixXd Xsig_aug = AugmentedSigmaPoints(x_in, P_in, std_a, std_yawdd);
  MatrixXd Xsig_pred = PredictSigmaPoints(Xsig_aug, delta_t);
  PredictMeanAndCovariance(x, P, Xsig_pred, weights);

  *Xsig_pred_out = Xsig_pred;

}

MatrixXd UKF::AugmentedSigmaPoints(VectorXd x,
                                  MatrixXd P,
                                  const double std_a, 
                                  const double std_yawdd)
{
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_ +1);

  //create augmented mean state
  x_aug.head(5) = x;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P;
  P_aug(5,5) = std_a*std_a;
  P_aug(6,6) = std_yawdd*std_yawdd;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  const double lambda = 3 - n_aug_;
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda+n_aug_) * L.col(i);
  }
  
  return Xsig_aug;

}

MatrixXd UKF::PredictSigmaPoints(MatrixXd Xsig_aug,
                                 const double delta_t)
{
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_  + 1);

  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }


  return Xsig_pred;
}

void UKF::PredictMeanAndCovariance(
  VectorXd* x_out,
  MatrixXd* P_out,
  MatrixXd Xsig_pred,
  VectorXd weights)
{
  VectorXd x = VectorXd(n_x_);
  MatrixXd P = MatrixXd(n_x_, n_x_);

  //predicted state mean
  x.fill(0.0);
  for (int i=0; i < 2*n_aug_ +1; i++)
  {
    x = x + weights(i) * Xsig_pred.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for(int i=0; i< 2* n_aug_ +1; i++)
  {
    VectorXd x_diff = Xsig_pred.col(i) -x;
    //angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2.f * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.f * M_PI;

    P = P + weights(i) * x_diff * x_diff.transpose();
  }

  *P_out = P;
  *x_out = x;

}
void UKF::PredictLidarMeasurement(MatrixXd* S_out,
                                  MatrixXd* Zsig_out,
                                  VectorXd* z_pred_out,
                                  MatrixXd Xsig_pred,
                                  VectorXd weights,
                                  const double std_laspx,
                                  const double std_laspy)
{

  const int dim_z = 2;
  MatrixXd Zsig = MatrixXd(dim_z, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(dim_z);
  MatrixXd S = MatrixXd(dim_z, dim_z);

  for(int i = 0; i< 2*n_aug_ + 1; i++)
  {
    double px = Xsig_pred(0, i);
    double py = Xsig_pred(1, i);
    Zsig(0, i) = px;
    Zsig(1, i) = py;

  }

  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights(i) * Zsig.col(i);
  }

  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(dim_z,dim_z);
  R <<    std_laspx*std_laspx, 0,
          0, std_laspy*std_laspy;
  S = S + R;

  *S_out = S;
  *Zsig_out = Zsig;
  *z_pred_out = z_pred;




}
void UKF::PredictRadarMeasurement(
                                  MatrixXd* S_out,
                                  MatrixXd* Zsig_out,
                                  VectorXd* z_pred_out,
                                  MatrixXd Xsig_pred,
                                  VectorXd weights,
                                  const double std_radr,
                                  const double std_radphi,
                                  const double std_radrd)
{
 

  const int dim_z = 3;
  MatrixXd Zsig = MatrixXd(dim_z, 2 * n_aug_ + 1);

  VectorXd z_pred = VectorXd(dim_z);

  MatrixXd S = MatrixXd(dim_z, dim_z);

  for (int i =0; i < 2*n_aug_ +1; i++){
      double p_x = Xsig_pred(0,i);
      double p_y = Xsig_pred(1,i);
      double v = Xsig_pred(2,i);
      double yaw = Xsig_pred(3,i);
      
      double v1 = cos(yaw)*v;
      double v2 = sin(yaw)*v;
      
      //measurement model
      Zsig(0,i) = sqrt(p_x*p_x  + p_y * p_y); //r
      Zsig(1,i) = atan2(p_y, p_x);              // phi
      Zsig(2,i) = (p_x*v1 + p_y*v2)/sqrt(p_x*p_x  + p_y * p_y); // r dot
      
  }

  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights(i) * Zsig.col(i);
  }

  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(dim_z,dim_z);
  R <<    std_radr*std_radr, 0, 0,
          0, std_radphi*std_radphi, 0,
          0, 0,std_radrd*std_radrd;
  S = S + R;

  *S_out = S;
  *Zsig_out = Zsig;
  *z_pred_out = z_pred;
}

void UKF::UpdateStateRadar(VectorXd* x_out, 
                          MatrixXd* P_out,
                          VectorXd z,
                          MatrixXd Zsig,
                          VectorXd z_pred,
                          MatrixXd Xsig_pred,
                          VectorXd x,
                          MatrixXd S,
                          MatrixXd P,
                          VectorXd weights
                          )
{
  //set measurement dimension, radar can measure r, phi, and r_dot
  int dim_z = 3;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, dim_z);
  Tc.fill(0.0);

  for(int i=0; i<2*n_aug_ + 1; i++)
  {
    VectorXd z_diff = (Zsig.col(i) - z_pred);
    while(z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while(z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    VectorXd x_diff = (Xsig_pred.col(i) - x);
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights(i) * x_diff*z_diff.transpose();
  }
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;
  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x = x + K * z_diff;
  P = P - K*S*K.transpose();

  *x_out = x;
  *P_out = P;
}

void UKF::UpdateStateLidar(VectorXd* x_out, 
                          MatrixXd* P_out,
                          VectorXd z,
                          MatrixXd Zsig,
                          VectorXd z_pred,
                          MatrixXd Xsig_pred,
                          VectorXd x,
                          MatrixXd S,
                          MatrixXd P,
                          VectorXd weights)
{
  //set measurement dimension, radar can measure r, phi, and r_dot
  int dim_z = 2;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, dim_z);
  Tc.fill(0.0);

  for(int i=0; i<2*n_aug_ + 1; i++)
  {
    VectorXd z_diff = (Zsig.col(i) - z_pred);
    VectorXd x_diff = (Xsig_pred.col(i) - x);
   
    Tc = Tc + weights(i) * x_diff*z_diff.transpose();
  }
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
  //update state mean and covariance matrix
  x = x + K * z_diff;
  P = P - K*S*K.transpose();

  *x_out = x;
  *P_out = P;
}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(
                      VectorXd* x,
                      MatrixXd* P,
                      MeasurementPackage meas_package,
                      MatrixXd Xsig_pred,
                      VectorXd weights,
                      const double std_laspx,
                      const double std_laspy) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  const int dim_z = 2; // px, py
  VectorXd z = meas_package.raw_measurements_;
  VectorXd x_in = *x;
  MatrixXd Zsig_buffer = MatrixXd(dim_z, 2 * n_aug_ + 1);
  MatrixXd S_buffer = MatrixXd(dim_z, dim_z);
  VectorXd z_pred_buffer = VectorXd(dim_z);
  MatrixXd P_in = *P;

  PredictLidarMeasurement(&S_buffer,
                          &Zsig_buffer,
                          &z_pred_buffer,
                          Xsig_pred,
                          weights,
                          std_laspx,
                          std_laspy);


  UpdateStateLidar(x, P, z, Zsig_buffer, z_pred_buffer, Xsig_pred, x_in, S_buffer,P_in, weights);

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(
                      VectorXd* x,
                      MatrixXd* P,
                      MeasurementPackage meas_package,
                      MatrixXd Xsig_pred,
                      VectorXd weights,
                      const double std_radr,
                      const double std_radphi,
                      const double stdradrd) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  const int dim_z = 3; // rho, phi, rho_dot
  VectorXd z = meas_package.raw_measurements_;
  VectorXd x_in = *x;

  MatrixXd Zsig_buffer = MatrixXd(dim_z, 2 * n_aug_ + 1);
  MatrixXd S_buffer = MatrixXd(dim_z, dim_z);
  VectorXd z_pred_buffer = VectorXd(dim_z);
  MatrixXd P_in = *P;
  PredictRadarMeasurement(
                          &S_buffer,
                          &Zsig_buffer,
                          &z_pred_buffer,
                          Xsig_pred,
                          weights,
                          std_radr,
                          std_radphi,
                          stdradrd);
  UpdateStateRadar(x, P, z, Zsig_buffer, z_pred_buffer, Xsig_pred, x_in, S_buffer,P_in, weights);


}

VectorXd UKF::ComputeWeights()
{
  VectorXd weights = VectorXd(2 * n_aug_ + 1);
  //set weights 
  weights[0] = 1.f/(lambda_ + n_aug_);
  for(int i=1; i < 2*n_aug_ +1; i++)
  {
    weights[i] = 0.5f/(lambda_ + n_aug_);
  }

  return weights;
}