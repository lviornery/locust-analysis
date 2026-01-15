close all
clear

method = ["butter","movav","polyfit","tvdiff"];
method = method(4);

butterOrder = 6;
nyquistFrac = 0.2;

mmn1 = 5;
mmn2 = 9;

polyorder = 9;

a1 = 0.1;
a2 = 0.1;

load("params_data.mat")

p = struct();
p = params_static(p);
p.dt = 1/2000;
p.tfinal = 15e-4; % seconds

exportTitles = {'Animal','Jump',...
    'Knee Angle (deg)','Spring Extension (m)',...
    'Knee Torque (N-m)','Spring Force (N)'};
exportmat = [];

for indivIdx = 1:7
    noffsets = sum(fileIndices(:,1)==indivIdx);
    offset_number = 0;
    offset_index = [];
    ext_m_k_full = [];
    force_m_k_full = [];
    phiRange = [inf,0];
    for idx = 1:length(fileIndices)
        if indivIdx
            if fileIndices(idx,1) ~= indivIdx
                continue
            end
        end
        offset_number = offset_number + 1;
        fileName = strcat("fulldata/a",num2str(fileIndices(idx,1)),"j",num2str(fileIndices(idx,2)),".csv");
        [t,vq1,vq2,vx,vy,vtheta,szr,mintthresh,maxtthresh,p] = system_extractJumpData(fileName,idx,p);
        
        z0 = de_bodyICs([vq1(1);vq2(1);vtheta(1);0;0;0],p);
        startX = z0(3);
        startY = z0(4);
        x = (vx - vx(1))*szr+startX;
        y = (vy - vy(1))*szr+startY;
        
        zeroOffset = [vx(1) - startX/szr,vy(1) - startY/szr];
        
        %% differentiation
        [Diff,Diff2] = matrixDiff(length(t),p.dt);
        
        if method == "butter"
            [b,a] = butter(butterOrder,nyquistFrac);
            vq1 = filtfilt(b,a,vq1);
            vq2 = filtfilt(b,a,vq2);
            x = filtfilt(b,a,x);
            y = filtfilt(b,a,y);
            vtheta = filtfilt(b,a,vtheta);
        
            vdq1 = Diff*vq1;
            vdq2 = Diff*vq2;
            dx = Diff*x;
            dy = Diff*y;
            vdtheta = Diff*vtheta;
            
            vddq1 = Diff2*vq1;
            vddq2 = Diff2*vq2;
            ddx = Diff2*x;
            ddy = Diff2*y;
            vddtheta = Diff2*vtheta;
        elseif method == "movav"    
            vdq1 = movmean(Diff*vq1,mmn1);
            vdq2 = movmean(Diff*vq2,mmn1);
            dx = movmean(Diff*x,mmn1);
            dy = movmean(Diff*y,mmn1);
            vdtheta = movmean(Diff*vtheta,mmn1);
            
            vddq1 = movmean(Diff2*vq1,mmn2);
            vddq2 = movmean(Diff2*vq2,mmn2);
            ddx = movmean(Diff2*x,mmn2);
            ddy = movmean(Diff2*y,mmn2);
            vddtheta = movmean(Diff2*vtheta,mmn2);
        elseif method == "polyfit"
            pq1 = polyfit(t,vq1,polyorder);
            pq2 = polyfit(t,vq2,polyorder);
            px = polyfit(t,x,polyorder);
            py = polyfit(t,y,polyorder);
            ptheta = polyfit(t,vtheta,polyorder);
            
            vdq1 = polyval(polyder(pq1),t);
            vdq2 = polyval(polyder(pq2),t);
            dx = polyval(polyder(px),t);
            dy = polyval(polyder(py),t);
            vdtheta = polyval(polyder(ptheta),t);
            
            vddq1 = polyval(polyder(polyder(pq1)),t);
            vddq2 = polyval(polyder(polyder(pq2)),t);
            ddx = polyval(polyder(polyder(px)),t);
            ddy = polyval(polyder(polyder(py)),t);
            vddtheta = polyval(polyder(polyder(ptheta)),t);
        elseif method == "tvdiff"
            [vdq1,vddq1] = regDiff(vq1,a1,a2,p.dt,[],1e10,1e-6,1e-10,5e-8);
            [vdq2,vddq2] = regDiff(vq2,a1,a2,p.dt,[],1e10,1e-6,1e-10,5e-8);
            [dx,ddx] = regDiff(x,a1,a2,p.dt,[],1e10,1e-6,1e-10,5e-8);
            [dy,ddy] = regDiff(y,a1,a2,p.dt,[],1e10,1e-6,1e-10,5e-8);
            [vdtheta,vddtheta] = regDiff(vtheta,a1,a2,p.dt,[],1e10,1e-6,1e-10,5e-8);
        end
        
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
            force_m_k(i) = -torque_m_k(i)/springArm;
            ext_m_k(i) = p.ExtensorDistance*((cos(phi(i) - p.ExtensorOffset))-(cos(phi(end) - p.ExtensorOffset)));
            force_s_k(i) = force_m_k(i)/ext_m_k(i);
        end

        maxextindex = find(vq2<p.MinKneeSpringAngle,1,"first");
        ext_m_k_full = [ext_m_k_full;ext_m_k(maxextindex:end)];
        force_m_k_full = [force_m_k_full;1*force_m_k(maxextindex:end)];
        new_offset_index = zeros(length(ext_m_k(maxextindex:end)),noffsets);
        new_offset_index(:,offset_number) = ones(size(ext_m_k(maxextindex:end)));
        offset_index = [offset_index;new_offset_index];

        if min(phi(maxextindex:end)) < phiRange(1)
            phiRange(1) = min(phi(maxextindex:end));
        end
        if max(phi(maxextindex:end)) > phiRange(2)
            phiRange(2) = max(phi(maxextindex:end));
        end

        n = length(t);
        exportAnimal = fileIndices(idx,1)*ones(n,1);
        exportJump = fileIndices(idx,2)*ones(n,1);
        exportAngle = rad2deg(phi(:));
        exportExtension = ext_m_k(:);
        exportTorque = -torque_m_k(:);
        exportForce = force_m_k(:);
        exportmat = [exportmat;exportAnimal,exportJump,exportAngle,exportExtension,exportTorque,exportForce];
    end
end

exportmat = [exportTitles;num2cell(exportmat)];
writecell(exportmat,'fig_6_export.csv')