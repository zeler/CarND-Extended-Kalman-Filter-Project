#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if(estimations.size() == 0 || estimations.size() != ground_truth.size()) {
      std::cout << "error" << std::endl;
  } else {

  // accumulate squared residuals
  for (int i=0; i < estimations.size(); ++i) {
    // ... your code here
    VectorXd v(4);
    
    v = estimations[i] - ground_truth[i];
    rmse = rmse.array() + (v.array() * v.array()); 
  }

  // calculate the mean
  rmse = rmse.array() / estimations.size();

  // calculate the squared root
  rmse = rmse.array().sqrt();

  }
  // return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // check division by zero
  if (px == 0 && py == 0) {
      std::cout << "Error, both zero" << std::endl;
  } else {
  // compute the Jacobian matrix
  
    float px2 = pow(px, 2);
    float py2 = pow(py, 2);
    float px2_plus_py2 = px2 + py2;
    float sqrt_px2_py2 = sqrt(px2_plus_py2);
    
    Hj << px / sqrt_px2_py2, py / sqrt_px2_py2, 0, 0,
          -py / px2_plus_py2, px / px2_plus_py2, 0, 0,
          py*(vx*py - vy*px) / pow(px2_plus_py2, 1.5f), px*(vy*px - vx*py) / pow(px2_plus_py2, 1.5f), px/sqrt_px2_py2, py / sqrt_px2_py2;
  }

  return Hj;
}
