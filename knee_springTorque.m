function springTorque = knee_springTorque(z,p)
%Takes a vector [q,dq]^T
%Returns the torque the spring exerts on the knee
q = z(1);
x = z(2);
dq = z(3);
dx = z(4);

%phi is the interior angle
phi = pi-q;
%calculate torque
springArm = p.ExtensorDistance*sin(phi - p.ExtensorOffset);
springForce = p.SpringConst*...
    p.ExtensorDistance*(cos(phi - p.ExtensorOffset)-p.springInitialExt);
springTorque = springArm*(springForce(z,p));