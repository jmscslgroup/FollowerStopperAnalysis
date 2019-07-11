% FollowerStopper Controller


clc
clf

 
%% SET PARAMETERS

w1 = 2;
w2 = 5.25;
w3 = 6.0;
time = data(:,1);
leader_xpos = data(:,2);
FollowerStopper_xpos = data(:,3);
leader_vel = data(:,4);
FollowerStopper_vel = data(:,5);

%relitive_xpos = leader_xpos - FollowerStopper_xpos;

%d_x  = relitive_xpos;
V_av = FollowerStopper_vel;
d_x = 10;
%V_av  = 22;
u    = 0.8;


d1 = (V_av.^2)./(2*u*(d_x-w1));
d2 = 1.0;
d3 = 0.5;
n = length(time);

dv = linspace(-10,0,n);
array = linspace(0,10,n);
w1_ = w1*ones(n,1);
w2_ = w2*ones(n,1);
w3_ = w3*ones(n,1);

 

 

%% CALCULATIONS
ksi1 = ones(n,n);
for i = 1:n
    ksi1(:,n) = w1 + 1/(2*d1(n)) * dv.^2;
end
ksi2 = w2 + 1/(2*d2) .* dv.^2;
ksi3 = w3 + 1/(2*d3) .* dv.^2;

 

 

%% GRAPH

plot(dv,ksi1,'b',dv,ksi2,'r',dv,ksi3,'g','LineWidth',2); hold on
plot(array,w1_,'b',array,w2_,'r',array,w3_,'g','LineWidth',2); hold off
axis([-5, 2, 0, 15]);
xlabel('\Delta v (m/s)'); ylabel('\Delta x (m)')
set(gca,'FontSize',18)
%lLine = animatedline('Color','k');%,'Marker','s','MarkerSize',10,'MarkerFaceColor','r');
%bLine = animatedline('Color','b');
%for n = 1:1:501
%    addpoints(lLine,data(n,1),data(n,2))
%    addpoints(bLine,ksi1)
%    drawnow
%    pause(0.05)
%end

