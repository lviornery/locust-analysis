function p = params_static(p)
%Returns a vector of system parameters to all the other happy little
%functions

p.lfemur = 22e-3;%.0508;%22e-3; %m - femur length
p.ltibia = 22e-3;%.0489;%22e-3; % m - tibia length
p.lbody = 30e-3; % m - body length - 10 mm to 45 mm
p.rbody = p.lbody*1/10; % m - body radius
p.hhip = 0;%-.02*p.lbody*(sqrt(2)/2); % m - height of hip relative to CoM
p.lhip = 0;%-.02*p.lbody*(sqrt(2)/2); % m - horizontal distance of hip relative to CoM
p.g = 9.81; %m/s^2
p.mbody = 1.1e-3; % kg - map
p.mfemur = p.mbody/20; %kg - map - NO
p.mtibia = p.mbody/100; %kg - map - NO
p.mtot = p.mbody + 2*p.mfemur + 2*p.mtibia;
p.Ibody = 1/12*p.mbody*(3*p.rbody^2+p.lbody^2); % kg-m - moment of inertia about major axis (apumes graphopper is a cylinder)
p.Ifemur = p.mfemur*p.lfemur^2/12;
p.Itibia = p.mtibia*p.ltibia^2/12;

p.SpringConst = 10;
p.springInitialExt = 0;
p.ExtensorDistance = 5.08e-3;%9.58e-4; %m
p.ExtensorOffset = deg2rad(0);%deg2rad(0); %rad

p.InitialHipAngle = deg2rad(120);
p.InitialLegAngle = deg2rad(180-11);
p.InitialBodyAngle = 0;

p.InitialHipSpeed = 0;
p.InitialLegSpeed = 0;
p.InitialBodySpeed = 0;

p.MaxHipAngle = deg2rad(180);
p.MinKneeAngle = deg2rad(180-90);
p.MinKneeSpringAngle = deg2rad(180-30);

% time parameters
p.dt = 1e-5; % seconds
p.tfinal = 15e-2; % seconds
tspan = 0:p.dt:p.tfinal;
p.kneeTorqueGrid = griddedInterpolant(tspan,-1*ones(size(tspan)));
p.hipTorqueGrid = griddedInterpolant(tspan,zeros(size(tspan)));
p.bodyXForceGrid = griddedInterpolant(tspan,zeros(size(tspan)));
p.bodyYForceGrid = griddedInterpolant(tspan,zeros(size(tspan)));
p.supportTorqueGrid = griddedInterpolant(tspan,zeros(size(tspan)));
end