% Chris Kreienkamp
% University of Arizona, CAT Vehicle
% July 3, 2019

%% stringStabilityTestData.m
% This file generate the lead vehicle data necessary to perform all of the
% string stability tests in the stringStabilityTest.slx Simulink file.

clear
clf
clc


%% OSCILLATORY
% Designed to collect transient data to understand how the follower behaves
% under non-constant headway and lead vehicle speed. Both 2.7 m/s and 4.5
% m/s speed fluctuations are tested. For the first half of the test the
% speed is fluctuated between 24.5 m/s and 21.9 m/s, with each speed being
% held for at least 30 seconds. For the second half of the test the speed
% is fluctuated between 24.5 m/s and 20.1 m/s with each speed being held
% for at least 30 seconds.
% Set FS desired velocity to 24.5 m/s.
osc_t = 30*9;
oscillatory(:,1) = linspace(0,osc_t,osc_t*100);
osc1 = 24.5*ones(osc_t*100/9,1);
osc2 = 21.9*ones(osc_t*100/9,1);
osc3 = 24.5*ones(osc_t*100/9,1);
osc4 = 21.9*ones(osc_t*100/9,1);
osc5 = 24.5*ones(osc_t*100/9,1);
osc6 = 20.1*ones(osc_t*100/9,1);
osc7 = 24.5*ones(osc_t*100/9,1);
osc8 = 20.1*ones(osc_t*100/9,1);
osc9 = 24.5*ones(osc_t*100/9,1);
oscillatory(:,2) = [osc1;osc2;osc3;osc4;osc5;osc6;osc7;osc8;osc9];


%% Low speed steps



%%