function u_out = fcn(max_speed,vel,max_accel,max_decel)
%#codegen
persistent y

if isempty(y)
    y=0;
end

% JMS: HACK timestep hacked at 0.05
dt=0.05;
if( y > max_speed + 1 )
    y = max(max_speed,y - abs(max_decel)*dt);
elseif( y < max_speed - 1)
    y = min(max_speed,y + max_accel*dt);
else
    y = double(max_speed);
end

if (y < 2 && max_speed > 2)
    y = 2;
elseif(y < 1 && max_speed > 1)
    y = 1;
end

u_out =min(max(y,vel-1.0),vel+2.0);