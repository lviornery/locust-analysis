function z0 = de_bodyICs(a0,p)
%Takes a vector [q1(0), q2(0), theta(0), dq1(0), dq2(0), dtheta(0)]^Te
%calculates x0, y0, dx0, and dy0
q1 = a0(1);
q2 = a0(2);
theta = a0(3);
dq1 = a0(4);
dq2 = a0(5);
dtheta = a0(6);

x0 = gen_xydxdyfromqs(dq1,dq2,dtheta,p.hhip,p.lfemur,p.lhip,p.ltibia,q1,q2,theta);
x = x0(1);
y = x0(2);
dx = x0(3);
dy = x0(4);
z0 = [q1;q2;x;y;theta;dq1;dq2;dx;dy;dtheta];
end

