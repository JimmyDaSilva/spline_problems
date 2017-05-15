#ifndef SPLINEPROBLEMS_EASYSCURVEPROFILE_HPP_
#define SPLINEPROBLEMS_EASYSCURVEPROFILE_HPP_

#include <ros/ros.h>
#include <std_msgs/Float32.h>
#include <eigen_conversions/eigen_kdl.h>
#include <spline_problems/matplotlibcpp.h>

class EasySCurveProfile{  
  // public functions
  public:
    // default constructor
    EasySCurveProfile();
    
    // constructor with initial parameters
    EasySCurveProfile(double s_init, double vi_init, double a_init, double s_final, double v_final, double a_final, double v_max, double a_max, double j_max);
    
    // change configuration
    void config (double s_init, double vi_init, double a_init, double s_final, double v_final, double a_final, double v_max, double a_max, double j_max);
    
    // config period
    void set_period ( double period );
    
    void compute_curves();
    void compute_next_step(double j);
    void compute_breaking();
  
    void plot_curves(std::string options ="");

    double compute_phase_distance(double time_in_phase, double j_phase, double phase_acc_start, double phase_vel_start);
    double compute_concave_distance(double phase_acc_start, double phase_vel_start, double phase_acc_final);
    double compute_concave_time(double phase_acc_start, double phase_acc_final);
    double compute_convexe_distance(double phase_acc_start, double phase_vel_start, double phase_acc_final);
    double compute_convexe_time(double phase_acc_start, double phase_acc_final);
    double compute_ramp_rise_distance(double phase_vel_start, double rise_vel_start, double rise_vel_final, double rise_acc_start, double rise_acc, double rise_acc_final);
    double compute_ramp_rise_time(double rise_vel_start, double rise_vel_final, double rise_acc_start, double rise_acc, double rise_acc_final );
    double compute_ramp_fall_distance(double phase_vel_start, double fall_vel_start, double fall_vel_final, double fall_acc_start, double fall_acc, double fall_acc_final);
    double compute_ramp_fall_time(double fall_vel_start, double fall_vel_final ,double fall_acc_start, double fall_acc, double fall_acc_final );
    double compute_cruise_distance(double cruise_vel, double phase_pos_start, double phase_pos_final);
    double compute_cruise_time(double cruise_vel, double phase_pos_start, double phase_pos_final);
  
  // public variables
  public:
    double period_;
    
  // private variables
  private:
    double v_max_;
    double a_max_;
    double j_max_;
    double si_;
    double vi_;
    double ai_;
    double sf_;
    double vf_;
    double af_;
    double break_dist_;
    double break_time_;
    
    std::vector<double> t_vect_, s_vect_, v_vect_, a_vect_, j_vect_;
};

#endif //SPLINEPROBLEMS_EASYSCURVEPROFILE_HPP_

