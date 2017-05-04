#include <ros/ros.h>
#include <std_msgs/Float32.h>
#include <eigen_conversions/eigen_kdl.h>
#include <spline_problems/matplotlibcpp.h>

double j_max = 5;
double v_max = 1;
double period = 0.001;
double qi = 0, qdi =0.5, qddi = 0.0;  
double qf = 2.0, qdf = 0.0, qddf = 0.0;

namespace plt = matplotlibcpp;

void curves_plot(std::vector<double>& t_vect, std::vector<double>& q_vect, std::vector<double>& v_vect, std::vector<double>& a_vect, std::vector<double>& j_vect, double t_rise);

int main(int argc, char** argv){
  ros::init(argc, argv, "tri_acc_profile");
  ros::NodeHandle nh;
  
  double t_rise = 2*std::sqrt((v_max-qdi)/j_max);
  
  std::vector<double> t_vect, q_vect, v_vect, a_vect, j_vect;
  
  // RISE CONCAVE
  for(double i = 0; i<=t_rise/2.0; i +=period)
    t_vect.push_back(i);
  
  int nb_rise_concave = t_vect.size();
  for(int i=0; i<nb_rise_concave;i++){
    j_vect.push_back(j_max);
    a_vect.push_back(j_max*t_vect[i]);
    v_vect.push_back(qdi+j_max / 2.0 * t_vect[i]* t_vect[i]);
    q_vect.push_back(qi+qdi*t_vect[i]+j_max/6.0*t_vect[i]*t_vect[i]*t_vect[i]);
  }
  
  // RISE CONVEXE
  double th = t_vect[t_vect.size()-1];
  double am = a_vect[t_vect.size()-1];
  double vh = v_vect[t_vect.size()-1];
  double sh = q_vect[t_vect.size()-1];
  for(double i = t_rise/2.0+period; i<=t_rise; i +=period)
    t_vect.push_back(i);
  
  int nb_rise_convexe = t_vect.size();
  
  double current_t;
  for(int i=nb_rise_concave; i<nb_rise_convexe;i++){
    current_t = t_vect[i] - th;
    j_vect.push_back(-j_max);
    a_vect.push_back(am - j_max*current_t);
    v_vect.push_back(vh+am*current_t - j_max*current_t*current_t/2.0);
    q_vect.push_back(sh+vh*current_t+am*current_t*current_t/2.0-j_max/6.0*current_t*current_t*current_t);
  }
  
  // CRUISE
  double ts = t_vect[t_vect.size()-1];
  double ss = q_vect[t_vect.size()-1];
  double t_fall = 2*std::sqrt((v_max-qdf)/j_max);;
  double breaking_distance = (v_max*v_max - qdf*qdf)*2/t_fall /j_max;
  double time_left = (qf - ss - breaking_distance)/v_max ;
  for(double i = t_rise+period; i<=t_rise+time_left; i +=period)
    t_vect.push_back(i);
  
  int nb_cruise = t_vect.size();
  
  for(int i=nb_rise_convexe; i<nb_cruise;i++){
    current_t = t_vect[i] - ts;
    j_vect.push_back(0.0);
    a_vect.push_back(0);
    v_vect.push_back(v_max);
    q_vect.push_back(ss+v_max*current_t);
  }
  
  // FALL CONVEXE
  for(double i = t_rise+time_left+period; i<=t_rise+time_left+t_fall/2.0; i +=period)
    t_vect.push_back(i);
  int nb_fall_convexe = t_vect.size();
  
  for(int i=nb_cruise; i<nb_fall_convexe;i++){
    current_t = t_vect[i] - t_rise - time_left;
    j_vect.push_back(-j_max);
    a_vect.push_back(- j_max*current_t);
    v_vect.push_back(v_max - j_max*current_t*current_t/2.0);
    q_vect.push_back(qf- breaking_distance+v_max*current_t-j_max/6.0*current_t*current_t*current_t);
  }
  
  // FALL CONCAVE
  am = a_vect[t_vect.size()-1];
  vh = v_vect[t_vect.size()-1];
  sh = q_vect[t_vect.size()-1];
  for(double i = t_rise+time_left+t_fall/2.0+period; i<=t_rise+time_left+t_fall; i +=period)
    t_vect.push_back(i);
  int nb_fall_concave = t_vect.size();
  
  for(int i=nb_fall_convexe; i<nb_fall_concave;i++){
    current_t = t_vect[i] - t_rise - time_left - t_fall/2.0;
    j_vect.push_back(j_max);
    a_vect.push_back(am + j_max*current_t);
    v_vect.push_back(vh + am*current_t + j_max*current_t*current_t/2.0);
    q_vect.push_back(sh+vh*current_t+ am/2.0*current_t*current_t+j_max/6.0*current_t*current_t*current_t);
  }

  // Plot
  curves_plot(t_vect, q_vect, v_vect, a_vect, j_vect, t_rise);

  ros::shutdown();
  return 1;
}

void curves_plot(std::vector<double>& t_vect, std::vector<double>& q_vect, std::vector<double>& v_vect, std::vector<double>& a_vect, std::vector<double>& j_vect, double t_rise){
  
  std::vector<double> line_t, line_val;
  
  // Plot position
  plt::subplot(2,2,1);
  plt::plot(t_vect, q_vect);
  plt::ylabel("position");
  line_t.clear();
  line_t.push_back(t_rise/2.0);
  line_t.push_back(t_rise/2.0);
  line_val.clear();
  line_val.push_back(0);
  line_val.push_back(qf);
  plt::plot(line_t, line_val,"r--");
  line_t.clear();
  line_t.push_back(t_rise);
  line_t.push_back(t_rise);
  line_val.clear();
  line_val.push_back(0);
  line_val.push_back(qf);
  plt::plot(line_t, line_val,"r--");
  
  // Plot velocity
  plt::subplot(2,2,2);
  plt::plot(t_vect, v_vect);
  plt::ylabel("velocity");
  line_t.clear();
  line_t.push_back(t_rise/2.0);
  line_t.push_back(t_rise/2.0);
  line_val.clear();
  line_val.push_back(0);
  line_val.push_back(v_max);
  plt::plot(line_t, line_val,"r--");
  line_t.clear();
  line_t.push_back(t_rise);
  line_t.push_back(t_rise);
  line_val.clear();
  line_val.push_back(0);
  line_val.push_back(v_max);
  plt::plot(line_t, line_val,"r--");
  
  // Plot acceleration
  plt::subplot(2,2,3);
  plt::plot(t_vect, a_vect);
  plt::ylabel("acceleration");
  line_t.clear();
  line_t.push_back(t_rise/2.0);
  line_t.push_back(t_rise/2.0);
  line_val.clear();
  line_val.push_back(0);
  line_val.push_back(j_max*t_rise/2.0);
  plt::plot(line_t, line_val,"r--");
  line_t.clear();
  line_t.push_back(t_rise);
  line_t.push_back(t_rise);
  line_val.clear();
  line_val.push_back(0);
  line_val.push_back(j_max*t_rise/2.0);
  plt::plot(line_t, line_val,"r--");
  
  // Plot jerk
  plt::subplot(2,2,4);
  plt::plot(t_vect, j_vect);
  plt::ylabel("jerk");
  line_t.clear();
  line_t.push_back(t_rise/2.0);
  line_t.push_back(t_rise/2.0);
  line_val.clear();
  line_val.push_back(-j_max);
  line_val.push_back(j_max);
  plt::plot(line_t, line_val,"r--");
  line_t.clear();
  line_t.push_back(t_rise);
  line_t.push_back(t_rise);
  line_val.clear();
  line_val.push_back(-j_max);
  line_val.push_back(j_max);
  plt::plot(line_t, line_val,"r--");
  
  plt::show();
  
}