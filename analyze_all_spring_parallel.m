close all
clear

butterOrder = 3;
nyquistFrac = 0.2;

mmn1 = 5;
mmn2 = 9;

polyorder = 6;

a1 = 0.1;
a2 = 0.1;

method = ["butter","movav","polyfit","tvdiff"];
method = method(4);

load("params_data.mat")
load("palattes\customPalatte.mat")

p = struct();
p = params_static(p);
p.dt = 1/2000;
p.tfinal = 15e-4; % seconds

fig = figure('Position', [10 10 5*150 800])
tlo = tiledlayout(1,2,"TileSpacing","compact")
tlr = tiledlayout(tlo,7,1,"TileSpacing","tight")
tlr.Layout.Tile = 1;
tll = tiledlayout(tlo,7,1,"TileSpacing","tight")
tll.Layout.Tile = 2;
subplotcount = 1;

post_rsqrs = [];

for indivIdx = 1:7
    nexttile(tlr,subplotcount)
    colororder(palatte)
    hold on
    nexttile(tll,subplotcount)
    colororder(palatte)
    hold on

    noffsets = sum(fileIndices(:,1)==indivIdx);
    offset_number = 0;
    offset_index = [];
    ext_m_k_full = [];
    force_m_k_full = [];
    extra_offset_index = [];
    extra_ext_m_k_full = [];
    extra_force_m_k_full = [];
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
        force_m_k_full = [force_m_k_full;force_m_k(maxextindex:end)];
        new_offset_index = zeros(length(ext_m_k(maxextindex:end)),noffsets);
        new_offset_index(:,offset_number) = ones(size(ext_m_k(maxextindex:end)));
        offset_index = [offset_index;new_offset_index];

        extra_ext_m_k_full = [extra_ext_m_k_full;ext_m_k(1:maxextindex)];
        extra_force_m_k_full = [extra_force_m_k_full;force_m_k(1:maxextindex)];
        extra_new_offset_index = zeros(length(ext_m_k(1:maxextindex)),noffsets);
        extra_new_offset_index(:,offset_number) = ones(size(ext_m_k(1:maxextindex)));
        extra_offset_index = [extra_offset_index;extra_new_offset_index];

        if min(phi(maxextindex:end)) < phiRange(1)
            phiRange(1) = min(phi(maxextindex:end));
        end
        if max(phi(maxextindex:end)) > phiRange(2)
            phiRange(2) = max(phi(maxextindex:end));
        end
        nexttile(tlr,subplotcount)
        scatter(ext_m_k,force_m_k,'x')
        nexttile(tll,subplotcount)
        scatter(rad2deg(phi),-torque_m_k,'x')
    end
    coefs = [offset_index,ext_m_k_full]\force_m_k_full;
    slope_index = noffsets+1;
    %dof_deflator is 1/(n-m), where n is the number of points and m is the number of coefficients in the fit
    dof_deflator = 1/(sum(offset_index(:))-slope_index);
    
    
    [X,SI] = sort(ext_m_k_full);
    Y = force_m_k_full(SI);
    
    yfmat = zeros(length(Y),noffsets);
    esqr = 0;
    xosqr = 0;
    yosqr = 0;
    esqr_post_offset = 0;
    yosqr_post_offset = 0;
    ax = nexttile(tlr,subplotcount)
    set(ax,'ColorOrderIndex',1)
    for i = 1:noffsets
        yfmat(:,i) = coefs(i)+coefs(slope_index)*X;
        plot(X,yfmat(:,i),LineWidth=1.5)
        %sum of squared errors of each data point's y-data
        esqr = esqr +...
            sum((coefs(i)+coefs(slope_index)*ext_m_k_full(offset_index(:,i)==1) - force_m_k_full(offset_index(:,i)==1)).^2);
        %sum of squared differences of each series' x-data from its own mean
        xosqr = xosqr +...
            sum((ext_m_k_full(offset_index(:,i)==1) - mean(ext_m_k_full(offset_index(:,i)==1))).^2);
        %sum of squared differences of each series' y-data from its own mean
        yosqr = yosqr +...
            sum((force_m_k_full(offset_index(:,i)==1) - mean(force_m_k_full(offset_index(:,i)==1))).^2);
        esqr_post_offset = esqr_post_offset +...
            sum((coefs(i)+coefs(slope_index)*extra_ext_m_k_full(extra_offset_index(:,i)==1) - extra_force_m_k_full(extra_offset_index(:,i)==1)).^2);
        yosqr_post_offset = yosqr_post_offset +...
            sum((extra_force_m_k_full(extra_offset_index(:,i)==1) - mean(extra_force_m_k_full(extra_offset_index(:,i)==1))).^2);
    end
    slope = coefs(slope_index)
    SEslope = sqrt((esqr/xosqr)*dof_deflator)
    rsqr = 1 - esqr/yosqr
    rsqr_post_offset = 1 - esqr_post_offset/yosqr_post_offset
    post_rsqrs = [post_rsqrs;rsqr_post_offset]

    ax = nexttile(tll,subplotcount)
    set(ax,'ColorOrderIndex',1)
    phiVec = linspace(phiRange(1),phiRange(2),100);
    springArmVec = p.ExtensorDistance*sin(phiVec - p.ExtensorOffset);
    ext_m_k_vec = p.ExtensorDistance*(cos(phiVec - p.ExtensorOffset));
    for i = 1:noffsets
        force_m_k_vec = coefs(i)+coefs(slope_index)*ext_m_k_vec;
        torque_m_k_vec = force_m_k_vec.*springArmVec;
        plot(rad2deg(phiVec),torque_m_k_vec,LineWidth=1.5)
    end
    if indivIdx == 7
        nexttile(tlr,subplotcount)
        text(0.6250*0.5e-4,0.4,{strcat("k \approx ",num2str(slope/1e3,3),"\pm",num2str(SEslope/1e3,2)," kN/m"),strcat("R^2 = ",num2str(rsqr,3))})
        ylim([0,0.5])
        nexttile(tll,subplotcount)
        ylim([0,1.5e-4])
    else
        nexttile(tlr,subplotcount)
        text(0.5e-4,2.5,{strcat("k \approx ",num2str(slope/1e3,3),"\pm",num2str(SEslope/1e3,1)," kN/m"),strcat("R^2 = ",num2str(rsqr,3))})
        ylim([0,3])
        nexttile(tll,subplotcount)
        ylim([0,1.25e-3])
    end
    nexttile(tlr,subplotcount)
    %xlim([0,0.8e-3])
    ylabel(strcat("Individual ",num2str(indivIdx)))
    nexttile(tll,subplotcount)
    xlim([0,100])

    subplotcount = subplotcount + 1;
end
ax = axes(tlo, 'Visible', 'off');
hold(ax, 'on')
scatter(ax,NaN,NaN,'x')
set(ax,'ColorOrderIndex',1)
plot(ax,NaN,NaN)
leg = legend(["Estimated Force/Torque";"Regression"],FontSize=10,Orientation="horizontal")
leg.IconColumnWidth = 10;
leg.Layout.Tile = 'north';
nexttile(tlr,7)
xlabel("$\Delta x\mathrm{ \, (m)}$",Interpreter="latex",FontSize=10)
nexttile(tll,7)
xlabel("$\theta_k \, \mathrm{(deg)}$",Interpreter="latex",FontSize=10)
ylabel(tlr,'Force (N)',FontSize=10);
ylabel(tll,'Torque (N-m)',FontSize=1);
saveas(fig,"Figures/Figure_7.png")