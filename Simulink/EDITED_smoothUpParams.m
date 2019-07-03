function u_out = fcn(max_speed,vel,max_accel,max_decel)
% #codegen
persistent y

if isempty(y)
    y=vel;
end

dt=0.01;
if( y > max_speed + max_accel*dt)
    y = max(max_speed,y - abs(max_decel)*dt);
elseif( y < max_speed - abs(max_decel)*dt)
    y = min(max_speed,y + max_accel*dt);
else
    y = max_speed;
end

u_out = y;