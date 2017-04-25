#include <spline_problems/probleme1.hpp>
#include <std_msgs/Float32.h>

// Probleme 1
// Compute the joint trajectory from q(0) = 1 to q(2) = 4 with null initial
// and final velocities and accelerations.

bool traj_done_;
double ti_, tf_, qi_, dqi_, ddqi_, qf_, dqf_, ddqf_, tc_, qc_, dqc_, ddqc_;

double current_time_, q_out_, dq_out_, ddq_out_;
void compute_traj();


int main(int argc, char** argv){
  ros::init(argc, argv, "probleme1");
  ros::NodeHandle nh;
  
  ros::Publisher q_pub = nh.advertise<std_msgs::Float32>("/splines/position",1);
  ros::Publisher dq_pub = nh.advertise<std_msgs::Float32>("/splines/velocity",1);
  ros::Publisher ddq_pub = nh.advertise<std_msgs::Float32>("/splines/acceleration",1);
  
  traj_done_ = false;
  ti_ = 0.0;
  tf_ = 1.0;
  qi_ = 0;
  dqi_ = 0.0;
  ddqi_ = 0.0;
  qf_ = 3;
  dqf_ = 0.0;
  ddqf_ = 0.0;
  
  dqc_ = 7;
  tc_ = (qi_-qf_+ tf_*dqc_)/dqc_;
  ddqc_ = (dqc_*dqc_)/(qi_-qf_+dqc_*tf_);
  
  ros::Duration r(0.01);
  current_time_ = ti_;
  q_out_ = 0.0;
  dq_out_ = 0.0;
  ddq_out_ = 0.0;
  std::cout << "tc : " << tc_ <<std::endl;
  std::cout << "dqc : " << dqc_ <<std::endl;
  std::cout << "ddqc : " << ddqc_ <<std::endl;
  
  std_msgs::Float32 pos_out_msg, vel_out_msg, acc_out_msg;
  while(ros::ok() && !traj_done_){
    compute_traj();
    
    r.sleep();
    current_time_ += 0.01;
    
    if (current_time_ > tf_)
      traj_done_ = true;
    
    pos_out_msg.data = q_out_;
    q_pub.publish(pos_out_msg);
    vel_out_msg.data = dq_out_;
    dq_pub.publish(vel_out_msg);
    acc_out_msg.data = ddq_out_;
    ddq_pub.publish(acc_out_msg);
  }
  
  ros::shutdown();
  return 0;
}

void compute_traj(){
  if((current_time_ >= 0) &&(current_time_ <= tc_)){
    q_out_ = qi_ + 0.5 * ddqc_ * current_time_ * current_time_;
    dq_out_ = ddqc_*current_time_ ;
    ddq_out_ = ddqc_;
    return;
  }
  if ((current_time_ > tc_) &&(current_time_<=tf_-tc_)){
    q_out_ = qi_ + ddqc_* tc_ * (current_time_ - tc_/2) ;
//     q_out_ = qi_ + ddqc_* tc_ * tc_ /2 + dqc_*(current_time_ - tc_);
    dq_out_ = ddqc_*tc_/2.0;
    ddq_out_ = 0;
    return;
  }
  if((current_time_ < tf_) &&(current_time_>tf_-tc_)){
    q_out_ = qf_ - 0.5 * ddqc_* (tf_ - current_time_) * (tf_ - current_time_) ;
    dq_out_ = ddqc_*(tf_-current_time_);
    ddq_out_ = -ddqc_;
    return;
  }
}