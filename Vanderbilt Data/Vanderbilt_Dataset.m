% Chris Kreienkamp
% University of Arizona, CAT Vehicle
% July 1, 2019

% This script uses a dataset provided by Vanderbilt which gives position,
% velocity, and acceleration values for each car in the ring-road
% experiment. The data will be extracted in order to determine a synthetic
% lead vehicle velocity that will imitate a human driver.

clc
clear 


%% LOAD DATA
% Possible experiments A-H. Different instructions given for different
% experiments. CAT Vehicle intervened to modify the natural traffic flow in
% experiments F, G, and H.
% Instruction I: "safely follow the vehilce in front as if in rush hour
% traffic" (A, B, C, D, E)
% Instruction II: "drive by the same instructions as before, but in
% addition place an emphasis on closing the gap to the vehicle in front,
% whenever safety permits" (F, G, H)

% Experiment E chosen as the synthetic lead vehicle velocity data because
% it features the most oscillatory traffic without intervention by the CAT
% Vehicle to modify the flow of cars.

filename_E = 'expE.csv';
Table_E = readtable(filename_E);
E = table2array(Table_E(:,8:11));
disp(Table_E(1,:))
Number_of_Rows_In_Table_E = height(Table_E);
Number_of_Cars_In_Exp_E = 19;
rowsPerCar_E = Number_of_Rows_In_Table_E / Number_of_Cars_In_Exp_E;


%% ORGANIZE DATA
time = E(:,1);
dist = E(:,2);
vel = E(:,3);
acc = E(:,4);

dist_car_1 = [time(1:rowsPerCar_E) dist(1:rowsPerCar_E)];
vel_car_1 = [time(1:rowsPerCar_E) vel(1:rowsPerCar_E)];
acc_car_1 = [time(1:rowsPerCar_E) acc(1:rowsPerCar_E)];
ave_vel_car_1 = mean(vel_car_1(:,2),'all')


%% ORGANIZE DATA FOR ANOTHER CAR
% carNumber = ??;
% for index = 1:rowsPerCar_E
%     i = index + rowsPerCar_E*(carNumber-1);
%     dist_car_??(index) = dist(i);
%     vel_car_??(index) = vel(i);
%     acc_car_??(index) = acc(i);
%     
% end

%% UNKNOWN USE
%t = 0.03 * [0:12480]';
%wave.time = t;
%wave.signals.values = dist_car_2;
%wave.signals.dimensions =1;
