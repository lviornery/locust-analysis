function pgrid = plot_jointPoints(z,p)
%Takes a vector [q1, q2, theta, dq1, dq2, dtheta]^T and
%calculates all system points as vectors of: body, hip, knee, foot
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

pgrid = gen_pointpositions(p.hhip,p.lfemur,p.lhip,p.ltibia,q1,q2,theta,x,y);
end
