#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;



  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction( MatrixXd* Xsig_pred_out,
                      VectorXd* x,
                      MatrixXd* P,
                      VectorXd weights,
                      const double std_a, 
                      const double std_yawdd,
                      double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  
  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  

  VectorXd ComputeWeights();
  void PredictMeanAndCovariance(
                                VectorXd* x_out,
                                MatrixXd* P_out,
                                MatrixXd Xsig_pred,
                                VectorXd weights);
  MatrixXd AugmentedSigmaPoints(VectorXd x,
                                MatrixXd P,
                                const double std_a, 
                                const double std_yawdd);
  MatrixXd PredictSigmaPoints(MatrixXd Xsig_aug,
                                 const double delta_t);
  void PredictRadarMeasurement(
                                  MatrixXd* S_out,
                                  MatrixXd* Zsig_out,
                                  VectorXd* z_pred_out,
                                  MatrixXd Xsig_pred,
                                  VectorXd weights,
                                  const double std_radr,
                                  const double std_radphi,
                                  const double stdradrd);
  void PredictLidarMeasurement(MatrixXd* S_out,
                                  MatrixXd* Zsig_out,
                                  VectorXd* z_pred_out,
                                  MatrixXd Xsig_pred,
                                  VectorXd weights,
                                  const double std_laspx,
                                  const double std_laspy);
  void UpdateStateRadar(VectorXd* x_out, 
                          MatrixXd* P_out,
                          VectorXd z,
                          MatrixXd Zsig,
                          VectorXd z_pred,
                          MatrixXd Xsig_pred,
                          VectorXd x,
                          MatrixXd S,
                          MatrixXd P,
                          VectorXd weights
                          );
  void UpdateStateLidar(VectorXd* x_out, 
                          MatrixXd* P_out,
                          VectorXd z,
                          MatrixXd Zsig,
                          VectorXd z_pred,
                          MatrixXd Xsig_pred,
                          VectorXd x,
                          MatrixXd S,
                          MatrixXd P,
                          VectorXd weights
                          );
  void UpdateLidar(
                  VectorXd* x,
                  MatrixXd* P,
                  MeasurementPackage meas_package,
                  MatrixXd Xsig_pred,
                  VectorXd weights,
                  const double std_laspx,
                  const double std_laspy);

  void UpdateRadar(
                  VectorXd* x,
                  MatrixXd* P,
                  MeasurementPackage meas_package,
                  MatrixXd Xsig_pred,
                  VectorXd weights,
                  const double std_radr,
                  const double std_radphi,
                  const double stdradrd); 
};

#endif /* UKF_H */
