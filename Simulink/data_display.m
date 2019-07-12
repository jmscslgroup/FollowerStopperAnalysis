% Chris Kreienkamp
% 07/03/19
% University of Arizona - CAT Vehicle


%% data_display.m

% This graphs the xposition-time, velocity-time, acceleration-time, and 
% yposition-xposition of Simulink-generated vehicles.

clc
clf


%% ORGANIZE DATA
%numberOfCars = 10;

t = xpos{1}.Values.time;
lead_xpos = xpos{1}.Values.data;
lead_ypos = ypos{1}.Values.data;
% lead_vel = vel.signals(1).values;
%lead_acc = acc.signals(1).values;
% for i = 2:numberOfCars
%     follower_xpos(:,i-1) = xpos.signals(i).values;
%     follower_ypos(:,i-1) = ypos.signals(i).values;
%     follower_vel(:,i-1) = vel.signals(i).values;
%     %follower_acc(:,i-1) = acc.signals(i).values;
% end


%% PLOT DATA
% figure(1)
% plot(t,lead_xpos,'r-'); hold on
% plot(t,follower_xpos); hold off
% axis([min(t) max(t) min(follower_xpos(:,numberOfCars-1)) max(lead_xpos)])
% set(gca,'FontSize',18)
% xlabel('Time (s)'); ylabel('Longitudinal Position (m)')
% legend('Lead','Follower','Location','Best')
% 
% figure(2)
% plot(t,lead_vel,'r-'); hold on
% plot(t, follower_vel); hold off
% axis([min(t) max(t) -0.05 max(lead_vel)+1.5])
% set(gca,'FontSize',18)
% xlabel('Time (s)'); ylabel('Longitudinal Velocity (m/s)')
% legend('Lead','Follower 1','Follower 2','Follower 3','Follower 4',...
%     'Follower 5','Follower 6','Follower 7','Follower 8','Follower 9',...
%     'Location','Best')

% figure(3)
% plot(t_acc,lead_acc,'r-'); hold on
% plot(t_acc,follower_acc,'b-'); hold off
% set(gca,'FontSize',18)
% xlabel('Time (s)'); ylabel('Longitudinal Acceleration (m/s/s)')
% legend('Lead','Follower','Location','Best')


%% PLOT MOVING CARS
figure(1)
y1=min(lead_ypos)-1;
y2=max(lead_ypos)+1;
xlabel('x-position (m)'); ylabel('y-position (m)')
set(gca,'Fontsize',18)
lLine = animatedline('Color','r','Marker','s','MarkerSize',10,'MarkerFaceColor','r');
% f1Line = animatedline('Color','b','Marker','s','MarkerSize',10,'MarkerFaceColor','b');
% f2Line = animatedline('Color','b','Marker','s','MarkerSize',10,'MarkerFaceColor','b');
% f3Line = animatedline('Color','b','Marker','s','MarkerSize',10,'MarkerFaceColor','b');
% f4Line = animatedline('Color','b','Marker','s','MarkerSize',10,'MarkerFaceColor','b');
% f5Line = animatedline('Color','b','Marker','s','MarkerSize',10,'MarkerFaceColor','b');
% f6Line = animatedline('Color','b','Marker','s','MarkerSize',10,'MarkerFaceColor','b');
% f7Line = animatedline('Color','b','Marker','s','MarkerSize',10,'MarkerFaceColor','b');
% f8Line = animatedline('Color','b','Marker','s','MarkerSize',10,'MarkerFaceColor','b');
% f9Line = animatedline('Color','g','Marker','s','MarkerSize',10,'MarkerFaceColor','g');
for n = 1:1:length(lead_xpos)
    axis([min(lead_xpos)-1 max(lead_xpos)+1 y1 y2])
    axis square
    clearpoints(lLine);
%     clearpoints(lLine); clearpoints(f1Line);  clearpoints(f2Line);
%      clearpoints(f3Line);  clearpoints(f4Line);  clearpoints(f5Line);
%       clearpoints(f6Line);  clearpoints(f7Line);  clearpoints(f8Line);
%        clearpoints(f9Line);
    addpoints(lLine,lead_xpos(n),lead_ypos(n))
%     addpoints(f1Line,follower_xpos(n,1),follower_ypos(n,1))
%     addpoints(f2Line,follower_xpos(n,2),follower_ypos(n,2))
%     addpoints(f3Line,follower_xpos(n,3),follower_ypos(n,3))
%     addpoints(f4Line,follower_xpos(n,4),follower_ypos(n,4))
%     addpoints(f5Line,follower_xpos(n,5),follower_ypos(n,5))
%     addpoints(f6Line,follower_xpos(n,6),follower_ypos(n,6))
%     addpoints(f7Line,follower_xpos(n,7),follower_ypos(n,7))
%     addpoints(f8Line,follower_xpos(n,8),follower_ypos(n,8))
%     addpoints(f9Line,follower_xpos(n,9),follower_ypos(n,9))
pause(0.00001)
end