function dz = de_ground(t,z,p)
%Takes a vector [q1, q2, x,y, theta, dq1, dq2, dx, dy, dtheta]^T and the time
%calculates torques and coordinate accelerations
%then constructs the vector [dq1, dq2, dx, dy, dtheta, ddq1, ddq2, ddx, ddy, ddtheta]^T

%calculate applied torque and solution to system
systemsol = de_ground_systemsol(t,z,p);
dz = [z(6:10);systemsol(1:5)];