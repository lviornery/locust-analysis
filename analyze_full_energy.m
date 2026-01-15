close all
clear

load("params_data.mat")

i_jump = zeros(1,length(fileIndices));
e_jump = zeros(1,length(fileIndices));
k_jump = [];
bl_jump = zeros(1,length(fileIndices));
ll_jump = zeros(1,length(fileIndices));
for indivIdx = unique(fileIndices(:,1))'
    offsets = zeros(1,sum(fileIndices(:,1)==indivIdx));
    offset_number = 0;
    ext_m_k_full = [];
    force_m_k_full = [];
    offset_index = [];
    for idx = 1:length(fileIndices)
        if fileIndices(idx,1) ~= indivIdx
            continue
        end
        i_jump(idx) = fileIndices(idx,1);
        offset_number = offset_number + 1;
        p = struct();
        p = params_static(p);
        p.ExtensorOffset = 0;
        p.dt = 1/2000;
        p.tfinal = 15e-4; % seconds
        
        angcutoff = deg2rad(90);
        
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
        
        [t,vq1,vq2,vx,vy,vtheta,szr,mintthresh,maxtthresh,p] = system_extractJumpData(fileName,idx,angcutoff,p);
        
        bl_jump(idx) = p.lbody;
        ll_jump(idx) = p.lfemur;
        
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
        % ddba = 800*ones(size(ddba));
        
        ze = [vq1(end),vq2(end),x(end),y(end),vtheta(end),vdq1(end),vdq2(end),dx(end),dy(end),vdtheta(end)];
        e_jump(idx) = system_get_ke(ze,p);
        
        phi = pi - vq2;
        torque_m_h = zeros(length(t),1);
        torque_m_k = zeros(length(t),1);
        torque_m_b = zeros(length(t),1);
        force_m_k = zeros(length(t),1);
        ext_m_k = zeros(length(t),1);
        force_s_k = zeros(length(t),1);
        for i = 1:length(t)
            zi = [vq1(i);vq2(i);x(i);y(i);vtheta(i);vdq1(i);vdq2(i);dx(i);dy(i);vdtheta(i);vddq1(i);vddq2(i);ddx(i);ddy(i);vddtheta(i)];
            [torque_m_h(i),torque_m_k(i),torque_m_b(i)] = system_get_torque(zi,p);
            springArm = p.ExtensorDistance*sin(phi(i) - p.ExtensorOffset);
            force_m_k(i) = torque_m_k(i)/springArm;
            ext_m_k(i) = cos(phi(i) - p.ExtensorOffset);
            force_s_k(i) = force_m_k(i)/ext_m_k(i);
        end
    
        maxextthresh = max(ext_m_k)*0.95;
        maxextindex = find(ext_m_k>=maxextthresh,1,"last");
        ext_m_k_full = [ext_m_k_full;ext_m_k(maxextindex:end)];
        force_m_k_full = [force_m_k_full;-1*force_m_k(maxextindex:end)];
        new_offset_index = zeros(length(ext_m_k(maxextindex:end)),length(offsets));
        new_offset_index(:,offset_number) = ones(size(ext_m_k(maxextindex:end)));
        offset_index = [offset_index;new_offset_index];
    end
    coefs = [offset_index,ext_m_k_full]\force_m_k_full;
    slope_index = offset_number+1;
    %dof_deflator is 1/(n-m), where n is the number of points and m is the number of coefficients in the fit
    dof_deflator = 1/(sum(offset_index(:))-slope_index);

    [X,SI] = sort(ext_m_k_full);
    Y = force_m_k_full(SI);

    yfmat = zeros(length(Y),offset_number);
    esqr = 0;
    xosqr = 0;
    for i = 1:offset_number
        yfmat(:,i) = coefs(i)+coefs(slope_index)*X;
        plot(X,yfmat(:,i))
        esqr = esqr +...
            sum((coefs(i)+coefs(slope_index)*ext_m_k_full(offset_index(:,i)==1) - force_m_k_full(offset_index(:,i)==1)).^2);
        xosqr = xosqr +...
            sum((ext_m_k_full(offset_index(:,i)==1) - mean(ext_m_k_full(offset_index(:,i)==1))).^2);
    end
    k_jump = [k_jump;coefs(slope_index)*ones(offset_number,1)];
end
figure
hold on
for idx = unique(i_jump)
    plot(k_jump(i_jump==idx),e_jump(i_jump==idx),'o',LineWidth=1)
end
hold off
xlabel('Spring constant')
ylabel('Jump kinetic energy')
legend([...
    "Individual 1";...
    "Individual 2";...
    "Individual 3";...
    "Individual 4";...
    "Individual 6";...
    "Individual 7";...
    "Individual 8"])

figure
hold on
for idx = unique(i_jump)
    mass = bodyMasses(bodyMasses(:,1)==idx,2) +...
        2*femurMasses(femurMasses(:,1)==idx,2) +...
        2*tibiaMasses(tibiaMasses(:,1)==idx,2);
    mass = mass*ones(size(e_jump(i_jump==idx)));
    plot(mass,e_jump(i_jump==idx),'o',LineWidth=1)
end
hold off
xlabel('Mass')
ylabel('Jump kinetic energy')
legend([...
    "Individual 1";...
    "Individual 2";...
    "Individual 3";...
    "Individual 4";...
    "Individual 6";...
    "Individual 7";...
    "Individual 8"])

figure
hold on
for idx = unique(i_jump)
    mass = bodyMasses(bodyMasses(:,1)==idx,2) +...
        2*femurMasses(femurMasses(:,1)==idx,2) +...
        2*tibiaMasses(tibiaMasses(:,1)==idx,2);
    mass = mass*ones(size(k_jump(i_jump==idx)));
    plot(mass,k_jump(i_jump==idx),'o',LineWidth=1)
end
hold off
xlabel('Mass')
ylabel('Spring constant')
legend([...
    "Individual 1";...
    "Individual 2";...
    "Individual 3";...
    "Individual 4";...
    "Individual 6";...
    "Individual 7";...
    "Individual 8"])

figure
hold on
for idx = unique(i_jump)
    plot(ll_jump(i_jump==idx),e_jump(i_jump==idx),'o',LineWidth=1)
end
hold off
xlabel('Body length')
ylabel('Jump kinetic energy')
legend([...
    "Individual 1";...
    "Individual 2";...
    "Individual 3";...
    "Individual 4";...
    "Individual 6";...
    "Individual 7";...
    "Individual 8"])

figure
hold on
for idx = unique(i_jump)
    plot(ll_jump(i_jump==idx),k_jump(i_jump==idx),'o',LineWidth=1)
end
hold off
xlabel('Body length')
ylabel('Spring constant')
legend([...
    "Individual 1";...
    "Individual 2";...
    "Individual 3";...
    "Individual 4";...
    "Individual 6";...
    "Individual 7";...
    "Individual 8"])