#include <spline_problems/easy_s_curve.hpp>

EasySCurveProfile::EasySCurveProfile ( double s_init, double vi_init, double a_init, double s_final, double v_final, double a_final, double v_max, double a_max, double j_max ) {
  si_ = s_init;
  sf_ = s_final;
  vi_ = vi_init;
  vf_ = v_final;
  ai_ = a_init;
  af_ = a_final;
  v_max_ = v_max;
  a_max_ = a_max;
  j_max_ = j_max;
  period_ = 0.00001;
}

EasySCurveProfile::EasySCurveProfile() {
  EasySCurveProfile(0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 1.0, 2.0, 5.0);
}

void EasySCurveProfile::config (double s_init, double vi_init, double a_init, double s_final, double v_final, double a_final, double v_max, double a_max, double j_max) {
  si_ = s_init;
  sf_ = s_final;
  vi_ = vi_init;
  vf_ = v_final;
  ai_ = a_init;
  af_ = a_final;
  period_ = 0.001;
  
  // Max vel acc and jerk need to be positive
  v_max_ = std::abs(v_max);
  a_max_ = std::abs(a_max);
  j_max_ = std::abs(j_max);
  
  // vf cannot be over vmax or below -vmax
  if (vf_>0)
    vf_ = std::min(vf_, v_max_);
  else
    vf_ = std::max(vf_, -v_max_);
  
  // af cannot be over amax or below -amax
  if (af_>0)
    af_ = std::min(af_, a_max_);
  else
    af_ = std::max(af_, -a_max_);
}

void EasySCurveProfile::set_period ( double period ) {
  period_ = period;
}

void EasySCurveProfile::compute_curves(){
  // Reset vector with initial conditions
  j_vect_.clear();
  j_vect_.push_back(0.0);
  a_vect_.clear();
  a_vect_.push_back(ai_);
  v_vect_.clear();
  v_vect_.push_back(vi_);
  s_vect_.clear();
  s_vect_.push_back(si_);
  t_vect_.clear();
  t_vect_.push_back(0.0);
  
  double distance_left = sf_-si_;
  compute_breaking();
  std::cout << "breaking_distance init : " <<break_dist_ <<std::endl;
  
  bool started_breaking = false;
  bool too_fast_on_start = ((vi_ + ai_*ai_/(2*j_max_) > v_max_) && (ai_>0)) ||( vi_ >v_max_);
  bool on_cruise = (vi_ == v_max_) && (ai_ == 0);
  
  while((break_dist_ <= distance_left) && (break_time_>=0)){
    // If v=vmax and a=0, keep going
    if (on_cruise)
      compute_next_step(0);
    else{
      // if vi > vmax
      if(too_fast_on_start){
        if(started_breaking){
          // a almost 0
          if(a_vect_[a_vect_.size()-1]+j_max_*period_>0){
            compute_next_step(0);
            on_cruise = true;
          }
          else
            compute_next_step(j_max_);
        }else{
          if((v_vect_[v_vect_.size()-1] - a_vect_[a_vect_.size()-1]*a_vect_[a_vect_.size()-1]/(2*j_max_) <= v_max_ ) && (a_vect_[a_vect_.size()-1]<0) ){
            compute_next_step(j_max_);
            started_breaking = true;
          }
          else{
            if(a_vect_[a_vect_.size()-1] > -a_max_)
              compute_next_step(-j_max_);
            else
              compute_next_step(0);
          }
        }
      }
      // vi <vmax
      else{
        // breaking to reach vmax
        if(started_breaking){
          if (a_vect_[a_vect_.size()-1]-j_max_*period_ < 0){
            compute_next_step(0);
            on_cruise = true;
          }
          else
            compute_next_step(-j_max_);
        }
        else{
          if((v_vect_[v_vect_.size()-1] + a_vect_[a_vect_.size()-1]*a_vect_[a_vect_.size()-1]/(2*j_max_) >= v_max_) && (a_vect_[a_vect_.size()-1]>0) ){
            compute_next_step(-j_max_);
            started_breaking = true;
          }
          else{
            if (a_vect_[a_vect_.size()-1] < a_max_)
              compute_next_step(j_max_);
            else
              compute_next_step(0);
          }
        }
      }
    }
    compute_breaking();
    distance_left = sf_ - s_vect_[s_vect_.size()-1];
  }
  
  std::cout << "switching! Time is : " <<t_vect_[t_vect_.size()-1] <<std::endl;
  
  if(break_dist_>1.1 * distance_left){
    std::cerr << "Cannot reach goal... Breaking HARD!" << std::endl;
    af_ = 0;
    vf_ = 0;
  }
  
  compute_breaking();
  double t_final = t_vect_[t_vect_.size()-1] + break_time_;
  std::cout << "break time is : " <<break_time_ <<std::endl;
  std::cout << "break distance is : " <<break_dist_ <<std::endl;
  
  while(t_vect_[t_vect_.size()-1] <=t_final){
    if ((vf_ >= v_vect_[v_vect_.size()-1] + (af_*af_-a_vect_[a_vect_.size()-1]*a_vect_[a_vect_.size()-1])/(2*j_max_)) && (a_vect_[a_vect_.size()-1]<0))
      compute_next_step(j_max_);
    else{
      if (a_vect_[a_vect_.size()-1] > -a_max_)
        compute_next_step(-j_max_);
      else
        compute_next_step(0);
    }
  }

}

void EasySCurveProfile::compute_breaking(){
  double ac = a_vect_[a_vect_.size()-1];
  double vc = v_vect_[v_vect_.size()-1];
  
  double t_ramp_fall = compute_ramp_fall_time(vc, vf_, ac, -a_max_, af_);
  double t_convexe_fall = compute_convexe_time(ac, -a_max_);
  double t_concave_fall = compute_concave_time(-a_max_,af_);
  
  if (t_ramp_fall <0){
    t_ramp_fall = 0;
    
    double vi_bis, vf_bis;
    if(af_>=0)
      vf_bis = vf_ - af_*af_/(2*j_max_);
    else
      vf_bis = vf_ + af_*af_/(2*j_max_);
    
    if(ac>=0)
      vi_bis = vc + ac*ac/(2*j_max_);
    else
      vi_bis = vc - ac*ac/(2*j_max_);
    
    double ideal_t = std::sqrt((vi_bis-vf_bis)/j_max_);
    t_concave_fall = ideal_t + af_/j_max_;
    t_convexe_fall = ideal_t + ac/j_max_;
  }
  
  if((t_concave_fall<0) || (t_concave_fall<0)){
    break_time_ = -1;
    break_dist_ = -1;
    return;
  }
  
  break_time_ = t_convexe_fall+t_ramp_fall+t_concave_fall;
  
  if (t_ramp_fall == 0)
    break_dist_ = compute_convexe_distance(ac,vc,ac-j_max_*t_convexe_fall) +compute_concave_distance(ac-j_max_*t_convexe_fall,vc+ac*t_convexe_fall-j_max_/2*t_convexe_fall*t_convexe_fall,af_);
  else
    break_dist_ = compute_convexe_distance(ac, vc, -a_max_) + compute_ramp_fall_distance(vc+ac*t_convexe_fall-j_max_/2*t_convexe_fall*t_convexe_fall,vc,vf_,ac,-a_max_,af_) + compute_concave_distance(-a_max_,vc+ac*t_convexe_fall-j_max_/2*t_convexe_fall*t_convexe_fall-t_ramp_fall*a_max_, af_) ;
  
}

void EasySCurveProfile::compute_next_step(double j){
  double last_a = a_vect_[a_vect_.size()-1];
  double last_v = v_vect_[v_vect_.size()-1];
  double last_s = s_vect_[s_vect_.size()-1];
  double last_t = t_vect_[t_vect_.size()-1];
  
  j_vect_.push_back(j);
  a_vect_.push_back(last_a + j*period_);
  v_vect_.push_back(last_v + last_a*period_ + j/2.0 * period_*period_);
  s_vect_.push_back(last_s + last_v*period_ + last_a/2.0*period_*period_ + j/6.0*period_*period_*period_);
  t_vect_.push_back(last_t + period_);
}

double EasySCurveProfile::compute_phase_distance(double time_in_phase, double j_phase, double phase_acc_start, double phase_vel_start){
  return phase_vel_start*time_in_phase+phase_acc_start/2.0*time_in_phase*time_in_phase+j_phase/6.0*time_in_phase*time_in_phase*time_in_phase;
}

double EasySCurveProfile::compute_concave_distance(double phase_acc_start, double phase_vel_start, double phase_acc_final){
  return compute_phase_distance(compute_concave_time(phase_acc_start, phase_acc_final), j_max_, phase_acc_start, phase_vel_start);
}

double EasySCurveProfile::compute_concave_time(double phase_acc_start, double phase_acc_final){
  return (phase_acc_final-phase_acc_start)/j_max_;
}

double EasySCurveProfile::compute_convexe_distance(double phase_acc_start, double phase_vel_start, double phase_acc_final){
  return compute_phase_distance(compute_convexe_time(phase_acc_start, phase_acc_final), -j_max_, phase_acc_start, phase_vel_start);
}

double EasySCurveProfile::compute_convexe_time(double phase_acc_start, double phase_acc_final){
  return (phase_acc_start-phase_acc_final)/j_max_;
}

double EasySCurveProfile::compute_ramp_rise_distance(double phase_vel_start, double rise_vel_start, double rise_vel_final, double rise_acc_start, double rise_acc, double rise_acc_final){
  return compute_phase_distance(compute_ramp_rise_time(rise_vel_start, rise_vel_final, rise_acc_start, rise_acc, rise_acc_final), 0 , rise_acc, phase_vel_start);
}

double EasySCurveProfile::compute_ramp_rise_time(double rise_vel_start, double rise_vel_final, double rise_acc_start, double rise_acc, double rise_acc_final ){
  return ((rise_vel_final-rise_vel_start)-(rise_acc*rise_acc-rise_acc_start*rise_acc_start/2-rise_acc_final*rise_acc_final/2.0)/j_max_)/rise_acc;
}

double EasySCurveProfile::compute_ramp_fall_distance(double phase_vel_start, double fall_vel_start, double fall_vel_final, double fall_acc_start, double fall_acc, double fall_acc_final){
  return compute_phase_distance(compute_ramp_fall_time(fall_vel_start, fall_vel_final, fall_acc_start, fall_acc, fall_acc_final), 0 , fall_acc, phase_vel_start);
}

double EasySCurveProfile::compute_ramp_fall_time(double fall_vel_start, double fall_vel_final ,double fall_acc_start, double fall_acc, double fall_acc_final ){
  return ((fall_vel_start-fall_vel_final)-(fall_acc*fall_acc-fall_acc_start*fall_acc_start/2-fall_acc_final*fall_acc_final/2)/j_max_)/(-fall_acc);
}

double EasySCurveProfile::compute_cruise_distance(double cruise_vel, double phase_pos_start, double phase_pos_final){
  return compute_phase_distance(compute_cruise_time(cruise_vel, phase_pos_start, phase_pos_final), 0, 0, cruise_vel);
}

double EasySCurveProfile::compute_cruise_time(double cruise_vel, double phase_pos_start, double phase_pos_final){
  return (phase_pos_final-phase_pos_start)/cruise_vel;
}

void EasySCurveProfile::plot_curves(std::string options) {
  if (options == ""){
    // Plot position
    matplotlibcpp::subplot(2,2,1);
    matplotlibcpp::plot(t_vect_, s_vect_);
    matplotlibcpp::ylabel("position");
    
    // Plot velocity
    matplotlibcpp::subplot(2,2,2);
    matplotlibcpp::plot(t_vect_, v_vect_);
    matplotlibcpp::ylabel("velocity");
    
    // Plot acceleration
    matplotlibcpp::subplot(2,2,3);
    matplotlibcpp::plot(t_vect_, a_vect_);
    matplotlibcpp::ylabel("acceleration");
    
    // Plot jerk
    matplotlibcpp::subplot(2,2,4);
    matplotlibcpp::plot(t_vect_, j_vect_);
    matplotlibcpp::ylabel("jerk");
  }
  else{
    // Plot position
    matplotlibcpp::subplot(2,2,1);
    matplotlibcpp::plot(t_vect_, s_vect_,options);
    matplotlibcpp::ylabel("position");
    
    // Plot velocity
    matplotlibcpp::subplot(2,2,2);
    matplotlibcpp::plot(t_vect_, v_vect_,options);
    matplotlibcpp::ylabel("velocity");
    
    // Plot acceleration
    matplotlibcpp::subplot(2,2,3);
    matplotlibcpp::plot(t_vect_, a_vect_,options);
    matplotlibcpp::ylabel("acceleration");
    
    // Plot jerk
    matplotlibcpp::subplot(2,2,4);
    matplotlibcpp::plot(t_vect_, j_vect_,options);
    matplotlibcpp::ylabel("jerk");
  }
}

int main(int argc, char** argv){
  ros::init(argc, argv, "easy_s_curve");
  ros::NodeHandle nh_param("~");
  double si, vi, ai, sf, vf, af, vmax, amax, jmax;
  nh_param.param<double>("s_init", si, 0);
  nh_param.param<double>("v_init", vi, 0);
  nh_param.param<double>("a_init", ai, 0);
  nh_param.param<double>("s_final", sf, 1);
  nh_param.param<double>("v_final", vf, 0);
  nh_param.param<double>("a_final", af, 0);
  nh_param.param<double>("v_max", vmax, 1);
  nh_param.param<double>("a_max", amax, 2);
  nh_param.param<double>("j_max", jmax, 10);

  EasySCurveProfile s_curve(si, vi, ai, sf, vf, af, vmax, amax, jmax);
  s_curve.compute_curves();
  s_curve.plot_curves();

  EasySCurveProfile s_curve2(si, vi, ai, sf, vf, af, vmax, amax, 1000);
  s_curve2.compute_curves();
  s_curve2.plot_curves("--r");
  
  matplotlibcpp::show();
  
  ros::shutdown();
  return 1;
}
