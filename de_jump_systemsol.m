function systemsol = de_jump_systemsol(t,z,p)
%Takes a vector [q1, q2, x, y, theta, dq1, dq2, dx, dy, dtheta]^T and the time
%calculates torques and coordinate accelerations
%then constructs the vector [dq1, dq2, dx, dy, dtheta, ddq1, ddq2, ddx, ddy, ddtheta]^T
q1 = z(1);
q2 = z(2);
x = z(3);
y = z(4);
theta = z(5);
dq1 = z(6);
dq2 = z(7);
dx = z(8);
dy = z(9);
dtheta = z(10);

%calculate applied torque and solution to system
kneeTorque = p.kneeTorqueGrid(t); % N-m
hipTorque = p.hipTorqueGrid(t); % N-m
bodyXForce = p.bodyXForceGrid(t); % N
bodyYForce = p.bodyYForceGrid(t); % N
supportTorque = p.supportTorqueGrid(t); % N-m
rhs = gen_jump_secondderivativesrhs(bodyXForce,bodyYForce,dq1,dq2,dtheta,dx,dy,p.g,p.hhip,hipTorque,kneeTorque,p.lfemur,p.lhip,p.ltibia,p.mbody,p.mfemur,p.mtibia,q1,q2,supportTorque,theta);
invm = gen_jump_secondderivativesinvm(p.Ibody,p.Ifemur,p.Itibia,p.hhip,p.lfemur,p.lhip,p.ltibia,p.mbody,p.mfemur,p.mtibia,q1,q2,theta,x,y);
systemsol = invm\rhs;