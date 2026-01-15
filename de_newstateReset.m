function [newState,newIC] = de_newstateReset(t,z,state,p)

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

newIC = z;

if abs(q1 - p.MaxHipAngle) < 1e-5
    newState = 3;
    newIC = [z(1:5);gen_hip_transition_sep(p.Ibody,p.Ifemur,p.Itibia,dq1,dq2,dtheta,dx,dy,p.hhip,p.lfemur,p.lhip,p.ltibia,p.mbody,p.mfemur,p.mtibia,q1,q2,theta)];
elseif abs(q2 - p.MinKneeAngle) < 1e-5
    newState = 3;
    newIC = [z(1:5);gen_knee_transition_sep(p.Ibody,p.Ifemur,p.Itibia,dq1,dq2,dtheta,dx,dy,p.hhip,p.lfemur,p.lhip,p.ltibia,p.mbody,p.mfemur,p.mtibia,q1,q2,theta)];
elseif (state==1)
    newState = 2;
elseif (state==2)
    systemsol = de_jump_systemsol(t,z,p);
    forces = systemsol(6:7);
    if abs(forces(2)) < 1e-5
        newState = 3;
    else
        newState = 2;
    end
else
    newState = 3;
end
