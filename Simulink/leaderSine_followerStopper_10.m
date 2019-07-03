% Chris Kreienkamp
% 06/27/19
% University of Arizona - CAT Vehicle


%% leaderSine_followerStopper_10.m
% This graphs 10 vehicles, with each vehicle starting 5 m behind the one
% immediately ahead. The lead vehicle travels at a sine wave velocity while
% the follower vehicles have their velocity controlled by followerStopper.

clear
clc
clf


%% ORGANIZE DATA
load leaderSine_followerStopper_10_xpos.mat
t_xpos  = leaderSine_followerStopper_10_xpos.time;
lead_xpos = leaderSine_followerStopper_10_xpos.signals(1).values;
follower1_xpos = leaderSine_followerStopper_10_xpos.signals(2).values;
follower2_xpos = leaderSine_followerStopper_10_xpos.signals(3).values;
follower3_xpos = leaderSine_followerStopper_10_xpos.signals(4).values;
follower4_xpos = leaderSine_followerStopper_10_xpos.signals(5).values;
follower5_xpos = leaderSine_followerStopper_10_xpos.signals(6).values;
follower6_xpos = leaderSine_followerStopper_10_xpos.signals(7).values;
follower7_xpos = leaderSine_followerStopper_10_xpos.signals(8).values;
follower8_xpos = leaderSine_followerStopper_10_xpos.signals(9).values;
follower9_xpos = leaderSine_followerStopper_10_xpos.signals(10).values;

load leaderSine_followerStopper_10_ypos.mat
t_ypos  = leaderSine_followerStopper_10_ypos.time;
lead_ypos = leaderSine_followerStopper_10_ypos.signals(1).values;
follower1_ypos = leaderSine_followerStopper_10_ypos.signals(2).values;
follower2_ypos = leaderSine_followerStopper_10_ypos.signals(3).values;
follower3_ypos = leaderSine_followerStopper_10_ypos.signals(4).values;
follower4_ypos = leaderSine_followerStopper_10_ypos.signals(5).values;
follower5_ypos = leaderSine_followerStopper_10_ypos.signals(6).values;
follower6_ypos = leaderSine_followerStopper_10_ypos.signals(7).values;
follower7_ypos = leaderSine_followerStopper_10_ypos.signals(8).values;
follower8_ypos = leaderSine_followerStopper_10_ypos.signals(9).values;
follower9_ypos = leaderSine_followerStopper_10_ypos.signals(10).values;

load leaderSine_followerStopper_10_vel.mat
t_vel  = leaderSine_followerStopper_10_vel.time;
lead_vel = leaderSine_followerStopper_10_vel.signals(1).values;
follower1_vel = leaderSine_followerStopper_10_vel.signals(2).values;
follower2_vel = leaderSine_followerStopper_10_vel.signals(3).values;
follower3_vel = leaderSine_followerStopper_10_vel.signals(4).values;
follower4_vel = leaderSine_followerStopper_10_vel.signals(5).values;
follower5_vel = leaderSine_followerStopper_10_vel.signals(6).values;
follower6_vel = leaderSine_followerStopper_10_vel.signals(7).values;
follower7_vel = leaderSine_followerStopper_10_vel.signals(8).values;
follower8_vel = leaderSine_followerStopper_10_vel.signals(9).values;
follower9_vel = leaderSine_followerStopper_10_vel.signals(10).values;

load leaderSine_followerStopper_10_acc.mat
t_acc  = leaderSine_followerStopper_10_acc.time;
lead_acc = leaderSine_followerStopper_10_acc.signals(1).values;
follower1_acc = leaderSine_followerStopper_10_acc.signals(2).values;
follower2_acc = leaderSine_followerStopper_10_acc.signals(3).values;
follower3_acc = leaderSine_followerStopper_10_acc.signals(4).values;
follower4_acc = leaderSine_followerStopper_10_acc.signals(5).values;
follower5_acc = leaderSine_followerStopper_10_acc.signals(6).values;
follower6_acc = leaderSine_followerStopper_10_acc.signals(7).values;
follower7_acc = leaderSine_followerStopper_10_acc.signals(8).values;
follower8_acc = leaderSine_followerStopper_10_acc.signals(9).values;
follower9_acc = leaderSine_followerStopper_10_acc.signals(10).values;


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
f4Line = animatedline('Color','b','Marker','s','MarkerSize',10,'MarkerFaceColor','b');
f5Line = animatedline('Color','b','Marker','s','MarkerSize',10,'MarkerFaceColor','b');
f6Line = animatedline('Color','b','Marker','s','MarkerSize',10,'MarkerFaceColor','b');
f7Line = animatedline('Color','b','Marker','s','MarkerSize',10,'MarkerFaceColor','b');
f8Line = animatedline('Color','b','Marker','s','MarkerSize',10,'MarkerFaceColor','b');
f9Line = animatedline('Color','g','Marker','s','MarkerSize',10,'MarkerFaceColor','g');
for n = 1:1:length(lead_xpos)
    axis([follower9_xpos(n)-20 lead_xpos(n)+20 y1 y2])
    clearpoints(lLine); clearpoints(f1Line);  clearpoints(f2Line);
     clearpoints(f3Line);  clearpoints(f4Line);  clearpoints(f5Line);
      clearpoints(f6Line);  clearpoints(f7Line);  clearpoints(f8Line);
       clearpoints(f9Line);
    addpoints(lLine,lead_xpos(n),lead_ypos(n))
    addpoints(f1Line,follower1_xpos(n),follower1_ypos(n))
    addpoints(f2Line,follower2_xpos(n),follower2_ypos(n))
    addpoints(f3Line,follower3_xpos(n),follower3_ypos(n))
    addpoints(f4Line,follower4_xpos(n),follower4_ypos(n))
    addpoints(f5Line,follower5_xpos(n),follower5_ypos(n))
    addpoints(f6Line,follower6_xpos(n),follower6_ypos(n))
    addpoints(f7Line,follower7_xpos(n),follower7_ypos(n))
    addpoints(f8Line,follower8_xpos(n),follower8_ypos(n))
    addpoints(f9Line,follower9_xpos(n),follower9_ypos(n))
    drawnow
end


%% PLOT DATA
figure(1)
plot(t_xpos,lead_xpos,'r-'); hold on
plot(t_xpos,follower1_xpos,t_xpos,follower2_xpos,t_xpos,follower3_xpos,...
    t_xpos,follower4_xpos,t_xpos,follower5_xpos,t_xpos,follower6_xpos,...
    t_xpos,follower7_xpos,t_xpos,follower8_xpos,t_xpos,follower9_xpos); hold off
set(gca,'FontSize',18)
xlabel('Time (s)'); ylabel('Longitudinal Position (m)')
legend('Lead','Location','Best')
axis([min(t_xpos) max(t_xpos) min(follower9_xpos) max(lead_xpos)])

figure(2)
plot(t_vel,lead_vel,'r-'); hold on
plot(t_vel,follower1_vel,t_xpos,follower2_vel,t_vel,follower3_vel,...
    t_vel,follower4_vel,t_vel,follower5_vel,t_vel,follower6_vel,...
    t_vel,follower7_vel,t_vel,follower8_vel,t_vel,follower9_vel); hold off
set(gca,'FontSize',18)
xlabel('Time (s)'); ylabel('Longitudinal Velocity (m/s)')
legend('Lead','Follower 1','Follower 2','Follower 3','Follower 4',...
    'Follower 5','Follower 6','Follower 7','Follower 8','Follower 9',...
    'Location','Best')
axis([min(t_vel) max(t_vel) min(follower9_vel) max(lead_vel)])

figure(3)
plot(t_acc,lead_acc,'r-'); hold on
plot(t_acc,follower1_acc,t_xpos,follower2_acc,t_acc,follower3_acc,...
    t_acc,follower4_acc,t_acc,follower5_acc,t_acc,follower6_acc,...
    t_acc,follower7_acc,t_acc,follower8_acc,t_acc,follower9_acc); hold off
set(gca,'FontSize',18)
xlabel('Time (s)'); ylabel('Longitudinal Acceleration (m/s/s)')
legend('Lead','Follower 1','Follower 2','Follower 3','Follower 4',...
    'Follower 5','Follower 6','Follower 7','Follower 8','Follower 9',...
    'Location','Best')
axis([min(t_acc) max(t_acc) -1 1])

