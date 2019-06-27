% Chris Kreienkamp
% 06/27/19
% University of Arizona - CAT Vehicle

% This graphs 2 vehicles, with the front vehicle starting at 10 m in front
% of the follower vehicle. The front vehicle, leadVehicle, travels at a
% constant velocity of 1 m/s while the follower vehicle, followerVehicle
% has its velocity controlled by followerStopper.

clear
clc
clf


%% ORGANIZE DATA
load leaderConstant_followerStopper_pos.mat
t_pos  = leaderConstant_followerStopper_pos.time;
lead_pos = leaderConstant_followerStopper_pos.signals(1).values;
follower_pos = leaderConstant_followerStopper_pos.signals(2).values;

load leaderConstant_followerStopper_vel.mat
t_vel  = leaderConstant_followerStopper_vel.time;
lead_vel = leaderConstant_followerStopper_vel.signals(1).values;
follower_vel = leaderConstant_followerStopper_vel.signals(2).values;

load leaderConstant_followerStopper_acc.mat
t_acc  = leaderConstant_followerStopper_acc.time;
lead_acc = leaderConstant_followerStopper_acc.signals(1).values;
follower_acc = leaderConstant_followerStopper_acc.signals(2).values;


%% PLOT DATA
figure(1)
plot(t_pos,lead_pos,'r-'); hold on
plot(t_pos,follower_pos,'b-'); hold off
set(gca,'FontSize',18)
xlabel('Time (s)'); ylabel('Longitudinal Position (m)')
legend('Lead','Follower','Location','Best')

figure(2)
plot(t_vel,lead_vel,'r-'); hold on
plot(t_vel,follower_vel,'b-'); hold off
set(gca,'FontSize',18)
xlabel('Time (s)'); ylabel('Longitudinal Velocity (m/s)')
legend('Lead','Follower','Location','Best')

figure(3)
plot(t_acc,lead_acc,'r-'); hold on
plot(t_acc,follower_acc,'b-'); hold off
set(gca,'FontSize',18)
xlabel('Time (s)'); ylabel('Longitudinal Acceleration (m/s)')
legend('Lead','Follower','Location','Best')
