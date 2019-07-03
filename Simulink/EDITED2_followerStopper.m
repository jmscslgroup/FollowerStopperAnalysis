function u_cmd = fcn(r,dx,dv,v_AV,...
    dx_min,dx_activate,decel)
%#codegen
% Safety controller, based on quadratic bands.
% Input: r = desired velocity (from other models)
%          dx = estimate of gap to vehicle ahead
%          dv = estimate of time derivative of dx (= velocity difference)
%         v_AV = velocity of AV
%      dx_min = minimum distance
% dx_activate = distance below which controller does something
%       decel = vector of three deceleration values for parabolas
% Out:  u_cmd = actual velocity commanded (always u_cmd<=U)
%
% Controller-specific parameters
dx_mid = (dx_min+dx_activate)/2; % mid distance, where v_lead is commanded
% Lead vehicle velocity
v_lead = v_AV+dv; % velocity of lead vehicle
v_lead = max(v_lead,0); % lead vehicle cannot go backwards
v = min(r,v_lead); % safety velocity cannot exceed desired velocity U
% Treatment of positive dv-values
dv = min(dv,0); % domains for dv>0 same as dv=0
% For given dv, dx-values of band boundaries
dx1 = dx_min+1/(2*decel(1))*dv.^2;
dx2 = dx_mid+1/(2*decel(2))*dv.^2;
dx3 = dx_activate+1/(2*decel(3))*dv.^2;

% Actual commanded velocity
dt = 0.01;
max_accel = 1.5;
smoothChange = v_AV + max_accel*dt;
if (dx <= dx1)
    u_cmd = 0;
elseif (dx <= dx2)
    u_cmd = min(v*(dx-dx1)/(dx2-dx1) , smoothChange);
elseif (dx <= dx3)
    u_cmd = min(v + (r-v)*(dx-dx2)/(dx3-dx2) , smoothChange);
else
    u_cmd = min(r , smoothChange);
end