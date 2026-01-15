function [tout,zout,ddzout,fout,stateout] = system_calculate(p)

%generate full-body ICs
a0 = params_IC(p);
z0 = de_bodyICs(a0,p);

%state = 1 means on the ground, state = 2 means jumping, state = 3
%means airborne
systemsol = de_jump_systemsol(0,z0,p);
force = systemsol(7);
if force > 0
    state = 2;
else
    state = 1;
end

if isequal(state,1)
    dz = de_ground(0,z0,p);
elseif isequal(state,2)
    dz = de_jump(0,z0,p);
else
    dz = de_air(0,z0,p);
end
dz = dz(6:10);

% Initialize output arrays
tout = 0;
zout = z0';
ddzout = dz';
stateout = state;
teout = [];
zeout = [];
ieout = [];
fout = force;
tstart = 0;

while tstart < p.tfinal
    % Initialize simulation time vector
    tspan = tstart:p.dt:p.tfinal+p.dt;
    % Tell ode45 what event function to use and set max step size to make sure
    % we don't miss a zero crossing - this goes in the loop so contactMode
    % can be updated
    options = odeset('Events', @(t,z) de_guard(t,z,state,p),'MaxStep',0.0001);
    
    % Simulate with ode45
    if isequal(state,1)
        [t,z,te,ze,ie] = ode45(@(t,z) de_ground(t,z,p),tspan,z0,options);
    elseif isequal(state,2)
        [t,z,te,ze,ie] = ode45(@(t,z) de_jump(t,z,p),tspan,z0,options);
    else
        [t,z,te,ze,ie] = ode45(@(t,z) de_air(t,z,p),tspan,z0,options);
    end
    
    % Sometimes the events function will record a nonterminal event if the
    % initial condition is a zero. We want to ignore this, so we will only
    % use the last row in the terminal state, time, and index.
    if ~isempty(ie)
        te = te(end,:);
        ze = ze(end,:);
        ie = ie(end,:);
    end
    
    % Log output
    nt = length(t);
    tout = [tout; t(2:nt)];
    zout = [zout; z(2:nt,:)];
    stateout = [stateout; state*ones(nt-1,1)];
    teout = [teout; te];
    zeout = [zeout; ze];
    ieout = [ieout; ie];
    for i = 2:nt
        if isequal(state,1)
            systemsol = de_ground_systemsol(t(i),z(i,:),p);
            f = systemsol(7);
            dz = systemsol(1:5);
        elseif isequal(state,2)
            systemsol = de_jump_systemsol(t(i),z(i,:),p);
            f = systemsol(7);
        else
            dz = de_air(t(i),z(i,:),p);
            f = 0;
        end
        dz = systemsol(1:5);
        ddzout = [ddzout; dz'];
        fout = [fout; f];
    end
    
    if isempty(ie)
        break; % abort if simulation has completed
    end
    
    [state,z0] = de_newstateReset(te,ze',state,p);
    tstart = t(end);
end

tvec = 0:p.dt:p.tfinal;
newZout = zeros(length(tvec),size(zout,2));
for i = 1:size(zout,2)
    newZout(:,i) = interp1(tout,zout(:,i),tvec);
end
zout = newZout;

newZout = zeros(length(tvec),size(ddzout,2));
for i = 1:size(ddzout,2)
    newZout(:,i) = interp1(tout,ddzout(:,i),tvec);
end
ddzout = newZout;

stateout = interp1(tout,stateout,tvec,"nearest");

fout = interp1(tout,fout,tvec);

tout = tvec';