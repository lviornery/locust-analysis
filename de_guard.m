function [constraintFcns,isterminal,direction] = de_guard(t,z,state,p)

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

constraintFcns = [q1 - p.MaxHipAngle;q2 - p.MinKneeAngle];

if(state==1)
    systemsol = de_ground_systemsol(t,z,p);
    forces = systemsol(6:8);
    constraintFcns = [constraintFcns;forces(2);forces(3)];
end

if(state==2)
    systemsol = de_jump_systemsol(t,z,p);
    forces = systemsol(6:7);
    constraintFcns = [constraintFcns;forces(2)];
end

isterminal = ones(length(constraintFcns),1);
direction = zeros(length(constraintFcns),1);
