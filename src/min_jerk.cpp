#include <ros/ros.h>
#include <std_msgs/Float32.h>
#include <eigen_conversions/eigen_kdl.h>
#include <spline_problems/matplotlibcpp.h>

double q_, qd_, qdd_, qddd_;

namespace plt = matplotlibcpp;

int main(int argc, char** argv){
  ros::init(argc, argv, "quintic_splines");
  ros::NodeHandle nh;
  
  double qi = 0.0, qdi = 20, qddi = 0.0;
  double tf = 0.5;
  double qf = 10, qdf = -10.0, qddf = 0.0;
  
  // A * coeffs = X
  Eigen::Matrix<double,6,6> A;
  A.setZero();
  A(0,0) = 1.0;
  A(1,1) = 1.0;
  A(2,2) = 2.0;
  A(3,0) = 1;
  A(3,1) = tf;
  A(3,2) = std::pow(tf,2);
  A(3,3) = std::pow(tf,3);
  A(3,4) = std::pow(tf,4);
  A(3,5) = std::pow(tf,5);
  A(4,1) = 1;
  A(4,2) = 2*tf;
  A(4,3) = 3*std::pow(tf,2);
  A(4,4) = 4*std::pow(tf,3);
  A(4,5) = 5*std::pow(tf,4);
  A(5,2) = 2;
  A(5,3) = 6*tf;
  A(5,4) = 12*std::pow(tf,2);
  A(5,5) = 20*std::pow(tf,3);
  
  Eigen::VectorXd vector_x, coeffs;
  vector_x.resize(6);
  coeffs.resize(6);
  vector_x(0) = qi;
  vector_x(1) = qdi;
  vector_x(2) = qddi;
  vector_x(3) = qf;
  vector_x(4) = qdf;
  vector_x(5) = qddf;
  
  coeffs = A.inverse() * vector_x;
  
  std::cout<< coeffs <<std::endl;
  
  std::vector<double> t_vect, q_vect, v_vect, a_vect, j_vect;
  for(double i = 0; i<=tf; i +=0.01)
    t_vect.push_back(i);
    
  for(int i=0; i<t_vect.size();i++){
    q_vect.push_back(coeffs(0) + coeffs(1)*t_vect[i] + coeffs(2)*std::pow(t_vect[i],2) + coeffs(3)*std::pow(t_vect[i],3) +coeffs(4)*std::pow(t_vect[i],4) +coeffs(5)*std::pow(t_vect[i],5)  );
    v_vect.push_back(coeffs(1) + 2*coeffs(2)*t_vect[i] + 3*coeffs(3)*std::pow(t_vect[i],2) + 4*coeffs(4)*std::pow(t_vect[i],3) +5*coeffs(5)*std::pow(t_vect[i],4));
    a_vect.push_back(2*coeffs(2) + 6*coeffs(3)*t_vect[i] + 12*coeffs(4)*std::pow(t_vect[i],2) + 20*coeffs(5)*std::pow(t_vect[i],3));
    j_vect.push_back(6*coeffs(3) + 24*coeffs(4)*t_vect[i] + 60*coeffs(5)*std::pow(t_vect[i],2));
  }
  
  plt::subplot(2,2,1);
  plt::plot(t_vect, q_vect);
  plt::ylabel("position");
  plt::subplot(2,2,2);
  plt::plot(t_vect, v_vect);
  plt::ylabel("velocity");
  plt::subplot(2,2,3);
  plt::plot(t_vect, a_vect);
  plt::ylabel("acceleration");
  plt::subplot(2,2,4);
  plt::plot(t_vect, j_vect);
  plt::ylabel("jerk");
  plt::show();
  

  ros::shutdown();
  return 1;
}