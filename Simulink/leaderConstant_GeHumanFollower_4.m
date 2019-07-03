% Chris Kreienkamp
% 06/27/19
% University of Arizona - CAT Vehicle


%% leaderConstant_GeHumanFollower_4.m

% This graphs 4 vehicles, with the lead vehicle starting at 20 m in front
% of the follower vehicle. The lead vehicle travels at a Constant velocity
% while the follower vehicles have their velocity controlled by the human
% follower model presented by Ge.

clear
clc
clf


%% ORGANIZE DATA
load leaderConstant_GeHumanFollower_4_xpos.mat
t_xpos  = leaderConstant_GeHumanFollower_4_xpos.time;
lead_xpos = leaderConstant_GeHumanFollower_4_xpos.signals(1).values;
follower1_xpos = leaderConstant_GeHumanFollower_4_xpos.signals(2).values;
follower2_xpos = leaderConstant_GeHumanFollower_4_xpos.signals(3).values;
follower3_xpos = leaderConstant_GeHumanFollower_4_xpos.signals(4).values;


load leaderConstant_GeHumanFollower_4_ypos.mat
t_ypos  = leaderConstant_GeHumanFollower_4_ypos.time;
lead_ypos = leaderConstant_GeHumanFollower_4_ypos.signals(1).values;
follower1_ypos = leaderConstant_GeHumanFollower_4_ypos.signals(2).values;
follower2_ypos = leaderConstant_GeHumanFollower_4_ypos.signals(3).values;
follower3_ypos = leaderConstant_GeHumanFollower_4_ypos.signals(4).values;


load leaderConstant_GeHumanFollower_4_vel.mat
t_vel  = leaderConstant_GeHumanFollower_4_vel.time;
lead_vel = leaderConstant_GeHumanFollower_4_vel.signals(1).values;
follower1_vel = leaderConstant_GeHumanFollower_4_vel.signals(2).values;
follower2_vel = leaderConstant_GeHumanFollower_4_vel.signals(3).values;
follower3_vel = leaderConstant_GeHumanFollower_4_vel.signals(4).values;


% load leaderConstant_GeHumanFollower_4_acc.mat
% t_acc  = leaderConstant_GeHumanFollower_4_acc.time;
% lead_acc = leaderConstant_GeHumanFollower_4_acc.signals(1).values;
% follower_acc = leaderConstant_GeHumanFollower_4_acc.signals(2).values;


% PLOT DATA
figure(1)
plot(t_xpos,lead_xpos,'r-'); hold on
plot(t_xpos,follower1_xpos,t_xpos,follower2_xpos,t_xpos,follower3_xpos); hold off
set(gca,'FontSize',18)
xlabel('Time (s)'); ylabel('Longitudinal Position (m)')
legend('Lead','Follower','Location','Best')

figure(2)
plot(t_vel,lead_vel,'r-'); hold on
plot(t_vel,follower1_vel,t_xpos,follower2_vel,t_vel,follower3_vel); hold off
set(gca,'FontSize',18)
xlabel('Time (s)'); ylabel('Longitudinal Velocity (m/s)')
legend('Lead','Follower 1','Follower 2','Follower 3','Location','Best')

% figure(3)
% plot(t_acc,lead_acc,'r-'); hold on
% plot(t_acc,follower_acc,'b-'); hold off
% set(gca,'FontSize',18)
% xlabel('Time (s)'); ylabel('Longitudinal Acceleration (m/s/s)')
% legend('Lead','Follower','Location','Best')


%% PLOT MOVING CARS
figure(4)
y1=min(lead_ypos)-1;
y2=max(lead_ypos)+1;
xlabel('x-position (m)'); ylabel('y-position (m)')
set(gca,'Fontsize',18)
lLine = animatedline('Color','r','Marker','s','MarkerSize',10,'MarkerFaceColor','r');
f1Line = animatedline('Color','b','Marker','s','MarkerSize',10,'MarkerFaceColor','b');
f2Line = animatedline('Color','b','Marker','s','MarkerSize',10,'MarkerFaceColor','b');
f3Line = animatedline('Color','b','Marker','s','MarkerSize',10,'MarkerFaceColor','b');
for n = 1:1:length(lead_xpos)
    axis([min(follower3_xpos)-0.5 max(lead_xpos)+0.5 y1 y2])
    clearpoints(lLine); clearpoints(f1Line);  clearpoints(f2Line);
     clearpoints(f3Line);  
    addpoints(lLine,lead_xpos(n),lead_ypos(n))
    addpoints(f1Line,follower1_xpos(n),follower1_ypos(n))
    addpoints(f2Line,follower2_xpos(n),follower2_ypos(n))
    addpoints(f3Line,follower3_xpos(n),follower3_ypos(n))
    drawnow
    %pause(0.01)
end