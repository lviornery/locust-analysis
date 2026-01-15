clear all
clc
close all

% save flag
savecommand = false;

%% variables
syms q1 dq1 ddq1 q2 dq2 ddq2 x dx ddx y dy ddy theta dtheta ddtheta kneeTorque hipTorque supportTorque bodyXForce bodyYForce lfemur ltibia hhip lhip...
    mbody mfemur mtibia Ibody Ifemur Itibia g PANIC real;
%phi = ground angle
%q1 = knee angle
%q2 = hip angle
%theta = jumper ground angle
%lf = length of femur
%lt = length of tibia
%hh = heightwise distance from hip to CoM
%lh = lengthwide distance from hip to CoM
%mo = mass of body
%ml = mass of leg segment
%Iz = moment of inertia of a body about its center
%g = gravity
psi = [q1;q2];
dpsi = [dq1;dq2];
ddpsi = [ddq1;ddq2];
z = [x;y;theta];
dz = [dx;dy;dtheta];
ddz = [ddx;ddy;ddtheta];
q = [psi;z];
dq = [dpsi;dz];
ddq = [ddpsi;ddz];
tau = [hipTorque;kneeTorque;bodyXForce;bodyYForce;supportTorque];
phi = -(theta+q1+q2);
dphi = -(dtheta+dq1+dq2);
mbodymat = [mbody,0,0,0,0,0;...
    0,mbody,0,0,0,0;...
    0,0,mbody,0,0,0;...
    0,0,0,PANIC,0,0;...
    0,0,0,0,PANIC,0;...
    0,0,0,0,0,Ibody];
mfemurmat = [mfemur,0,0,0,0,0;...
    0,mfemur,0,0,0,0;...
    0,0,mfemur,0,0,0;...
    0,0,0,PANIC,0,0;...
    0,0,0,0,PANIC,0;...
    0,0,0,0,0,Ifemur];
mtibiamat = [mtibia,0,0,0,0,0;...
    0,mtibia,0,0,0,0;...
    0,0,mtibia,0,0,0;...
    0,0,0,PANIC,0,0;...
    0,0,0,0,PANIC,0;...
    0,0,0,0,0,Itibia];
%% frames
gos = [1,0,0,lhip;...
    0,1,0,hhip;...
    0,0,1,0;...
    0,0,0,1];
gsh = [cos(q1),-sin(q1),0,0;...
    sin(q1),cos(q1),0,0;...
    0,0,1,0;...
    0,0,0,1];
ghfc = [1,0,0,lfemur/2;...
    0,1,0,0;...
    0,0,1,0;...
    0,0,0,1];
ghk = [cos(q2),-sin(q2),0,lfemur;...
    sin(q2),cos(q2),0,0;...
    0,0,1,0;...
    0,0,0,1];
gktc = [1,0,0,ltibia/2;...
    0,1,0,0;...
    0,0,1,0;...
    0,0,0,1];
gka = [1,0,0,ltibia;...
    0,1,0,0;...
    0,0,1,0;...
    0,0,0,1];
gwo = [cos(theta),-sin(theta),0,x;...
    sin(theta),cos(theta),0,y;...
    0,0,1,0;...
    0,0,0,1];
gwc = [1,0,0,0;...
    0,1,0,0;...
    0,0,1,0;...
    0,0,0,1];
gac = [cos(phi),-sin(phi),0,0;...
    sin(phi),cos(phi),0,0;...
    0,0,1,0;...
    0,0,0,1];
gofc = simplify(gos*gsh*ghfc,'Steps',10);
gotc = simplify(gos*gsh*ghk*gktc,'Steps',10);
%% constraint
%grasp side
B = [eye(2);zeros(4,2)];
gco = inv(gwc)*gwo;
Adgco = simplify(rigidbody_tform2adjoint(gco));
Gs = -Adgco'*B;
Jbop = [cos(theta),sin(theta),0;...
    -sin(theta),cos(theta),0;...
    zeros(3);...
    0,0,1];
gbart = simplify(Gs'*Jbop);

%hand Jacobian side
gsa = gsh*ghk*gka;
gsc = simplify(gsa*gac);
Adgscinv = simplify(rigidbody_tform2adjoint(inv(gsc)));
%s jacobian of the fingertips
Jssa = rigidbody_spatialjacobian(gsa,psi);
Jh = simplify(B'*Adgscinv*Jssa);

%constraint
A = simplify([-Jh, gbart],'Steps',10);
%acceleration constraint
dA = sym(zeros(size(A)));
for i = 1:length(A)
    dA(:,i) = simplify(jacobian(A(:,i),q)*dq,'Steps',10);
end
%% Lagrange... sorta
gwfc = simplify(gwo*gofc,'Steps',10);
gwtc = simplify(gwo*gos*gsh*ghk*gktc,'Steps',10);
pfc = gwfc*[0;0;0;1];
ptc = gwtc*[0;0;0;1];
po = gwo*[0;0;0;1];
PE = simplify(pfc(2)*mfemur*g+ptc(2)*mtibia*g + po(2)*mbody*g,'Steps',10);

%kinetic energy - get equivalent mass matrix first
Jbofc = simplify(rigidbody_bodyjacobian(gofc,q),'Steps',10);
Jbotc = simplify(rigidbody_bodyjacobian(gotc,q),'Steps',10);
Jbwo = simplify(rigidbody_bodyjacobian(gwo,q),'Steps',10);
Adgofcinv = simplify(rigidbody_tform2adjoint(inv(gofc)),'Steps',10);
Adgotcinv = simplify(rigidbody_tform2adjoint(inv(gotc)),'Steps',10);
topleft = simplify(Jbofc'*mfemurmat*Jbofc + Jbotc'*mtibiamat*Jbotc,'Steps',10);
topright = simplify(Jbofc'*mfemurmat*Adgofcinv*Jbwo + Jbotc'*mtibiamat*Adgotcinv*Jbwo,'Steps',10);
bottomleft = simplify(Jbwo'*Adgofcinv'*mfemurmat*Jbofc + Jbwo'*Adgotcinv'*mtibiamat*Jbotc,'Steps',10);
bottomright = simplify(Jbwo'*Adgofcinv'*mfemurmat*Adgofcinv*Jbwo + ...
    Jbwo'*Adgotcinv'*mtibiamat*Adgotcinv*Jbwo + Jbwo'*mbodymat*Jbwo,'Steps',10);
mbar = topleft+topright+bottomleft+bottomright;

nbar = simplify(jacobian(PE,q)','Steps',10);

invm = [mbar,A';A,zeros(size(A,1))];
%Lagrange's equation
systemrhs = [tau-nbar-rigidbody_coriolismatrix(mbar,q,dq)*dq;-dA*dq];

minirhs = rigidbody_coriolismatrix(mbar,q,dq)*dq + nbar;

psol = [po, gwo*gos*gsh*[0;0;0;1], gwo*gos*gsh*ghk*[0;0;0;1], gwo*gos*gsh*ghk*gka*[0;0;0;1]];

pbalt = gwc*inv(gos*gsh*ghk*gka*gac)*[0;0;0;1];
%% Simplify
invm = simplify(invm,'IgnoreAnalyticConstraints',true,'Steps',10);
systemrhs = simplify(systemrhs,'IgnoreAnalyticConstraints',true,'Steps',10);
xyequiv = simplify(pbalt(1:2),'IgnoreAnalyticConstraints',true,'Steps',10);
dxdyequiv = simplify(jacobian(xyequiv,[psi;theta])*[dpsi;dtheta],'IgnoreAnalyticConstraints',true,'Steps',10);
ddxddyequiv = simplify(jacobian(dxdyequiv,[psi;theta;dpsi;dtheta])*[dpsi;dtheta;ddpsi;ddtheta],'IgnoreAnalyticConstraints',true,'Steps',10);
psol = simplify(psol(1:2,:),'IgnoreAnalyticConstraints',true,'Steps',10);
KE = simplify(.5*dq'*mbar*dq,'IgnoreAnalyticConstraints',true,'Steps',10);
%% Invert to find torque
%dy = Ydq
%ddy = dYdq + Yddq
%dY = 0, ddy = Yddq
%approx_torque_tilde = mtilde*ddy+ctilde*dy + ntilde
%assuming only nonzero q1, q2, theta applied forces,
%approx_torque_tilde = [t_q1;t_q2;t_theta]
Y = [1,0,0,0,0;0,1,0,0,0;0,0,0,0,1];
H = [A;Y]\[zeros(2,3);eye(3)];
H = simplify(H,'IgnoreAnalyticConstraints',true,'Steps',10);
dH = sym(zeros(size(H)));
for i = 1:size(H,2)
    dH(:,i) = simplify(jacobian(H(:,i),q)*dq,'Steps',10);
end

%ddq method
dyr = Y*dq;
Hddyr = ddq - dH*dyr;
mtilde = H'*mbar;
ctilde = H'*mbar*dH+H'*rigidbody_coriolismatrix(mbar,q,dq)*H;
ntilde = H'*nbar;
approx_torque_ddq = mtilde*Hddyr+ctilde*dyr+ntilde;
approx_torque_ddq = simplify(approx_torque_ddq,'IgnoreAnalyticConstraints',true,'Steps',10);

%dq method
dyr = Y*dq;
ddyr = Y*ddq;
mtilde = H'*mbar*H;
ctilde = H'*mbar*dH+H'*rigidbody_coriolismatrix(mbar,q,dq)*H;
ntilde = H'*nbar;
approx_torque_dq = mtilde*ddyr+ctilde*dyr+ntilde;
approx_torque = simplify(0.5*(approx_torque_dq+approx_torque_ddq),'IgnoreAnalyticConstraints',true,'Steps',10);
%% Impulse functions
%hip joint end of RoM, no ground contact
Ai = [1,0,0,0,0];
AdagT = Ai*inv(mbar)*Ai';
AdagT = inv(mbar)*Ai'*inv(AdagT);
AdagT = simplify(AdagT,'IgnoreAnalyticConstraints',true,'Steps',10);
hip_transition_sep = dq - AdagT*Ai*dq;

%knee joint end of RoM, no ground contact
Ai = [0,1,0,0,0];
AdagT = Ai*inv(mbar)*Ai';
AdagT = inv(mbar)*Ai'*inv(AdagT);
AdagT = simplify(AdagT,'IgnoreAnalyticConstraints',true,'Steps',10);
knee_transition_sep = dq - AdagT*Ai*dq;
%% Create functions
matlabFunction(hip_transition_sep,'File','gen_hip_transition_sep')
matlabFunction(knee_transition_sep,'File','gen_knee_transition_sep')
matlabFunction(approx_torque,'File','gen_approx_torque')
matlabFunction(A,'File','gen_A')
matlabFunction(dA,'File','gen_dA')
matlabFunction(systemrhs,'File','gen_jump_secondderivativesrhs')
matlabFunction(invm,'File','gen_jump_secondderivativesinvm')
matlabFunction([xyequiv;dxdyequiv],'File','gen_xydxdyfromqs')
matlabFunction(ddxddyequiv,'File','gen_ddxddyfromqs')
matlabFunction(psol,'File','gen_pointpositions')
matlabFunction(KE,'File','gen_ke')
% matlabFunction(PE,'File','gen_pe')
%% Create constraints on ground
%y-velocity restricted to zero
A = [A;...
    0,0,0,1,0];
dA = sym(zeros(size(A)));
for i = 1:length(A)
    dA(:,i) = simplify(jacobian(A(:,i),q)*dq,'Steps',10);
end
%% Generate system equations on ground
invm = [mbar,A';A,zeros(size(A,1))];
%Lagrange's equation
systemrhs = [tau-nbar-rigidbody_coriolismatrix(mbar,q,dq)*dq;-dA*dq];
%% Simplify
invm = simplify(invm,'IgnoreAnalyticConstraints',true,'Steps',10);
systemrhs = simplify(systemrhs,'IgnoreAnalyticConstraints',true,'Steps',10);
%% Create functions
matlabFunction(systemrhs,'File','gen_ground_secondderivativesrhs')
matlabFunction(invm,'File','gen_ground_secondderivativesinvm')

clear all