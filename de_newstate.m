function newState = de_newstate(t,z,state,p)

if (state==1)
    newState = 2;
elseif (state ==2)
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
