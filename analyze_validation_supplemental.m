close all
clear

idx = 28;

load("params_data.mat")

p = struct();
p = params_static(p);
p.ExtensorOffset = deg2rad(5);
p.dt = 1/2000;
p.tfinal = 15e-4; % seconds

angcutoff = deg2rad(80);

doLineVideo = false;
doWireframe = false;
doOverlay = true;

a1q1 = 0.01;
a2q1 = 0.1;
a1q2 = 0.01;
a2q2 = 0.1;
a1x = 0.01;
a2x = 0.1;
a1y = 0.01;
a2y = 0.1;
a1theta = 0.01;
a2theta = 0.1;

fileName = strcat("fulldata/a",num2str(fileIndices(idx,1)),"j",num2str(fileIndices(idx,2)),".csv");
videoName = strcat("fulldata/Locust ",num2str(fileIndices(idx,1))," jump ",num2str(fileIndices(idx,2)),".mp4");

[t,vq1,vq2,vx,vy,vtheta,szr,mintthresh,maxtthresh,p] = system_extractJumpData(fileName,idx,angcutoff,p);

p.MinKneeAngle = min(vq2);

z0 = de_bodyICs([vq1(1);vq2(1);vtheta(1);0;0;0],p);

startX = z0(3);
startY = z0(4);
x = (vx - vx(1))*szr+startX;
y = (vy - vy(1))*szr+startY;

zeroOffset = [vx(1) - startX/szr,vy(1) - startY/szr];

%% differentiation
[vdq1,vddq1] = regDiff(vq1,a1q1,a2q1,p.dt,[],1e10,1e-6,1e-10,5e-8);
[vdq2,vddq2] = regDiff(vq2,a1q2,a2q2,p.dt,[],1e10,1e-6,1e-10,5e-8);
[dx,ddx] = regDiff(x,a1x,a2x,p.dt,[],1e10,1e-6,1e-10,5e-8);
[dy,ddy] = regDiff(y,a1y,a2y,p.dt,[],1e10,1e-6,1e-10,5e-8);
[vdtheta,vddtheta] = regDiff(vtheta,a1theta,a2theta,p.dt,[],1e10,1e-6,1e-10,5e-8);

torque_m_h = zeros(length(t),1);
torque_m_k = zeros(length(t),1);
torque_m_b = zeros(length(t),1);
force_m_k = zeros(length(t),1);
ext_m_k = zeros(length(t),1);
force_s_k = zeros(length(t),1);
predictForward = zeros(length(t),10);
cons_viol = zeros(length(t),4);
for i = 1:length(t)
    zi = [vq1(i);vq2(i);x(i);y(i);vtheta(i);vdq1(i);vdq2(i);dx(i);dy(i);vdtheta(i);vddq1(i);vddq2(i);ddx(i);ddy(i);vddtheta(i)];
    [torque_m_h(i),torque_m_k(i),torque_m_b(i)] = system_get_torque(zi,p);
    cons_viol(i,1:2) = gen_A(p.lfemur,p.ltibia,zi(1),zi(2),zi(5),zi(3),zi(4))*zi(6:10);
    cons_viol(i,3:4) = gen_dA(zi(6),zi(7),zi(10),zi(8),zi(9),p.lfemur,p.ltibia,zi(1),zi(2),zi(5))*zi(11:15);
    phi = pi-vq2(i);
    springArm = p.ExtensorDistance*sin(phi - p.ExtensorOffset);
    force_m_k(i) = torque_m_k(i)/springArm;
    ext_m_k(i) = cos(phi - p.ExtensorOffset);
    force_s_k(i) = force_m_k(i)/ext_m_k(i);
end

p.hipTorqueGrid = griddedInterpolant(t,torque_m_h);
p.kneeTorqueGrid = griddedInterpolant(t,torque_m_k);
p.supportTorqueGrid = griddedInterpolant(t,torque_m_b);
for i = 1:length(t)
    predictForward(i,:) = de_jump(t(i),[vq1(i),vq2(i),x(i),y(i),vtheta(i),vdq1(i),vdq2(i),dx(i),dy(i),vdtheta(i)],p);
end

figure
timeVX = axes;
hold on
figure
timeVY = axes;
hold on

ke = .5*p.mtot*(dx.^2+dy.^2);
% forceMag = sqrt(force(:,1).^2 + force(:,2).^2 );
keAll = ke(end);
% plot(t,ddx,t,ddy)
plot(timeVX,t,dx,t,cumtrapz(p.dt,ddx)+dx(1))
plot(timeVY,t,dy,t,cumtrapz(p.dt,ddy)+dy(1))

title(timeVX,"Integrated acceleration and speed for each trajectory in X")
title(timeVY,"Integrated acceleration and speed for each trajectory in Y")