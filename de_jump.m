function dz = de_jump(t,z,p)
%Takes a vector [q1, q2, x,y, theta, dq1, dq2, dx, dy, dtheta]^T and the time
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
systemsol = de_jump_systemsol(t,z,p);
dz = [[dq1;dq2;dx;dy;dtheta];systemsol(1:5)];