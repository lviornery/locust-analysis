function a0 = params_IC(p)
%initial conditions
q10 = p.InitialHipAngle;
q20 = p.InitialLegAngle;
theta0 = p.InitialBodyAngle;
beta = rad2deg(atan2(p.lfemur*sin(q10+theta0)-p.ltibia*sin(q10+q20+theta0),-p.lfemur*cos(q10+theta0)-p.ltibia*cos(q10+q20+theta0)))
dq10 = p.InitialHipSpeed;
dq20 = p.InitialLegSpeed;
dtheta0 = p.InitialBodySpeed;
a0 = [q10,q20,theta0,dq10,dq20,dtheta0];
end

