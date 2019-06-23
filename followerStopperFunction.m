function [v_cmd,ksi1,ksi2,ksi3] = followerStopperFunction(xRel,vRel,v,vDesired)

%% SET PARAMETERS
% minimum distance parameters
w1 = 4.5;
w3 = 6.0;
% maximum acceleration parameters
a1 = 1.5;
a3 = 0.5;
% others
w2 = mean([w1,w3],'all');
a2 = mean([a1,a3],'all');
vLead = v + vRel;
vStar = min(v,0);
nu = min(max(vLead,0),vDesired);


%% CALCULATIONS
% ksi values
ksi1 = w1 + 1/(2*a1) * vStar^2;
ksi2 = w2 + 1/(2*a2) * vStar^2;
ksi3 = w3 + 1/(2*a3) * vStar^2;
% command velocity
if xRel <= ksi1
    v_cmd = 0;
elseif xRel <= ksi2
    v_cmd = nu*(xRel-ksi1)/(ksi2-ksi1);
elseif xRel <= ksi3
    v_cmd = nu + (vDesired - nu)*(xRel-ksi2)/(ksi3-ksi2);
else
    v_cmd = vDesired;
end

end