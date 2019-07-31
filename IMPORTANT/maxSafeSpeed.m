% Chris Kreienkamp
% University of Arizona, CAT Vehicle
% July 29, 2019
% Calculation of the max safe speed of the AV

clear
clc


%% PARAMETERS
G = 9.80665;
psi = 1;
ksi = 81;
delta = 1.158;
a_max = 3.53;
a_dmax = -7.66;
k = -G/a_dmax;
v_lead = 0;
v_AV = linspace(0,35,10000);


%% CALCULATION
deltaVSStar = max(0, 1/2/k/a_dmax.*(v_lead^2-k.*v_AV.^2));
y = -ksi + psi + deltaVSStar + v_AV.*(1-a_max/a_dmax)*delta...
    + a_max/2*(1-a_max/a_dmax)*delta^2;
% y = -ksi + psi + deltaVSStar + v_AV.*delta;

%% PLOT
plot(v_AV,y)
set(gca,'FontSize',18)
xlabel('v_{AV}'); ylabel('y')
