% Chris Kreienkamp
% 06/27/19
% University of Arizona - CAT Vehicle


%% leaderRandom_followerStopper.m

% This graphs 2 vehicles, with the lead vehicle starting at 10 m in front
% of the follower vehicle. The lead vehicle travels at a random velocity
% while the follower vehicle has its velocity controlled by 
% followerStopper.

clear
clc
clf


%% ORGANIZE DATA
load leaderRandom_followerStopper_xpos.mat
t_xpos  = leaderRandom_followerStopper_xpos.time;
lead_xpos = leaderRandom_followerStopper_xpos.signals(1).values;
follower_xpos = leaderRandom_followerStopper_xpos.signals(2).values;

load leaderRandom_followerStopper_ypos.mat
t_ypos  = leaderRandom_followerStopper_ypos.time;
lead_ypos = leaderRandom_followerStopper_ypos.signals(1).values;
follower_ypos = leaderRandom_followerStopper_ypos.signals(2).values;

load leaderRandom_followerStopper_vel.mat
t_vel  = leaderRandom_followerStopper_vel.time;
lead_vel = leaderRandom_followerStopper_vel.signals(1).values;
follower_vel = leaderRandom_followerStopper_vel.signals(2).values;

load leaderRandom_followerStopper_acc.mat
t_acc  = leaderRandom_followerStopper_acc.time;
lead_acc = leaderRandom_followerStopper_acc.signals(1).values;
follower_acc = leaderRandom_followerStopper_acc.signals(2).values;


%% PLOT DATA
figure(1)
plot(t_xpos,lead_xpos,'r-'); hold on
plot(t_xpos,follower_xpos,'b-'); hold off
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
xlabel('Time (s)'); ylabel('Longitudinal Acceleration (m/s/s)')
legend('Lead','Follower','Location','Best')


%% PLOT MOVING CARS
figure(4)
y1=min(follower_ypos)-1;
y2=max(follower_ypos)+1;
xlabel('x-position (m)'); ylabel('y-position (m)')
set(gca,'Fontsize',18)
lLine = animatedline('Color','r','Marker','s','MarkerSize',10,'MarkerFaceColor','r');
fLine = animatedline('Color','b','Marker','s','MarkerSize',10,'MarkerFaceColor','b');
for n = 1:1:length(follower_xpos)
    axis([min(follower_xpos)-0.5 max(lead_xpos)+0.5 y1 y2])
    clearpoints(fLine)
    clearpoints(lLine)
    addpoints(lLine,lead_xpos(n),lead_ypos(n))
    addpoints(fLine,follower_xpos(n),follower_ypos(n))
    drawnow
    pause(0.25)
end
