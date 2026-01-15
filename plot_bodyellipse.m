function ellipsecoordinates = plot_bodyellipse(z,p)
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

t = linspace(0,2*pi,100);
a = p.lbody/2;
b = p.rbody;
ellipsecoordinates = zeros(2,length(t));
ellipsecoordinates(1,:) = x + a*cos(t)*cos(theta) - b*sin(t)*sin(theta);
ellipsecoordinates(2,:)  = y + b*sin(t)*cos(theta) + a*cos(t)*sin(theta);
end

