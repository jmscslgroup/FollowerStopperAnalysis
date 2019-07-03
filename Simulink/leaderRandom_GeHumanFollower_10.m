% Chris Kreienkamp
% 06/27/19
% University of Arizona - CAT Vehicle


%% leaderRandom_GeHumanFollower_4.m

% This graphs 4 vehicles, with the lead vehicle starting at 20 m in front
% of the follower vehicle. The lead vehicle travels at a random velocity
% while the follower vehicles have their velocity controlled by the human
% follower model presented by Ge.

clear
clc
clf


%% ORGANIZE DATA
load leaderRandom_GeHumanFollower_4_xpos.mat
t_xpos  = leaderRandom_GeHumanFollower_4_xpos.time;
lead_xpos = leaderRandom_GeHumanFollower_4_xpos.signals(1).values;
follower1_xpos = leaderRandom_GeHumanFollower_4_xpos.signals(2).values;
follower2_xpos = leaderRandom_GeHumanFollower_4_xpos.signals(3).values;
follower3_xpos = leaderRandom_GeHumanFollower_4_xpos.signals(4).values;
follower4_xpos = leaderRandom_GeHumanFollower_4_xpos.signals(5).values;
follower5_xpos = leaderRandom_GeHumanFollower_4_xpos.signals(6).values;
follower6_xpos = leaderRandom_GeHumanFollower_4_xpos.signals(7).values;
follower7_xpos = leaderRandom_GeHumanFollower_4_xpos.signals(8).values;
follower8_xpos = leaderRandom_GeHumanFollower_4_xpos.signals(9).values;
follower9_xpos = leaderRandom_GeHumanFollower_4_xpos.signals(10).values;

load leaderRandom_GeHumanFollower_4_ypos.mat
t_ypos  = leaderRandom_GeHumanFollower_4_ypos.time;
lead_ypos = leaderRandom_GeHumanFollower_4_ypos.signals(1).values;
follower1_ypos = leaderRandom_GeHumanFollower_4_ypos.signals(2).values;
follower2_ypos = leaderRandom_GeHumanFollower_4_ypos.signals(3).values;
follower3_ypos = leaderRandom_GeHumanFollower_4_ypos.signals(4).values;
follower4_ypos = leaderRandom_GeHumanFollower_4_ypos.signals(5).values;
follower5_ypos = leaderRandom_GeHumanFollower_4_ypos.signals(6).values;
follower6_ypos = leaderRandom_GeHumanFollower_4_ypos.signals(7).values;
follower7_ypos = leaderRandom_GeHumanFollower_4_ypos.signals(8).values;
follower8_ypos = leaderRandom_GeHumanFollower_4_ypos.signals(9).values;
follower9_ypos = leaderRandom_GeHumanFollower_4_ypos.signals(10).values;

load leaderRandom_GeHumanFollower_4_vel.mat
t_vel  = leaderRandom_GeHumanFollower_4_vel.time;
lead_vel = leaderRandom_GeHumanFollower_4_vel.signals(1).values;
follower1_vel = leaderRandom_GeHumanFollower_4_vel.signals(2).values;
follower2_vel = leaderRandom_GeHumanFollower_4_vel.signals(3).values;
follower3_vel = leaderRandom_GeHumanFollower_4_vel.signals(4).values;
follower4_vel = leaderRandom_GeHumanFollower_4_vel.signals(5).values;
follower5_vel = leaderRandom_GeHumanFollower_4_vel.signals(6).values;
follower6_vel = leaderRandom_GeHumanFollower_4_vel.signals(7).values;
follower7_vel = leaderRandom_GeHumanFollower_4_vel.signals(8).values;
follower8_vel = leaderRandom_GeHumanFollower_4_vel.signals(9).values;
follower9_vel = leaderRandom_GeHumanFollower_4_vel.signals(10).values;

% load leaderRandom_GeHumanFollower_4_acc.mat
% t_acc  = leaderRandom_GeHumanFollower_4_acc.time;
% lead_acc = leaderRandom_GeHumanFollower_4_acc.signals(1).values;
% follower_acc = leaderRandom_GeHumanFollower_4_acc.signals(2).values;


% PLOT DATA
figure(1)
plot(t_xpos,lead_xpos,'r-'); hold on
plot(t_xpos,follower1_xpos,t_xpos,follower2_xpos,t_xpos,follower3_xpos,...
    t_xpos,follower4_xpos,t_xpos,follower5_xpos,t_xpos,follower6_xpos,...
    t_xpos,follower7_xpos,t_xpos,follower8_xpos,t_xpos,follower9_xpos); hold off
set(gca,'FontSize',18)
xlabel('Time (s)'); ylabel('Longitudinal Position (m)')
legend('Lead','Follower','Location','Best')

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
f4Line = animatedline('Color','b','Marker','s','MarkerSize',10,'MarkerFaceColor','b');
f5Line = animatedline('Color','b','Marker','s','MarkerSize',10,'MarkerFaceColor','b');
f6Line = animatedline('Color','b','Marker','s','MarkerSize',10,'MarkerFaceColor','b');
f7Line = animatedline('Color','b','Marker','s','MarkerSize',10,'MarkerFaceColor','b');
f8Line = animatedline('Color','b','Marker','s','MarkerSize',10,'MarkerFaceColor','b');
f9Line = animatedline('Color','g','Marker','s','MarkerSize',10,'MarkerFaceColor','g');
for n = 1:1:length(lead_xpos)
    axis([min(follower9_xpos)-0.5 max(lead_xpos)+0.5 y1 y2])
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
    %pause(0.01)
end