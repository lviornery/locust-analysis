close all
clear

load("params_data.mat")
load("palattes\customPalatte.mat")

method = ["butter","movav","polyfit","tvdiff"];
method = method(4);

butterOrder = 3;
nyquistFrac = 0.2;

mmn1 = 5;
mmn2 = 9;

polyorder = 9;

a1 = 0.1;
a2 = 0.1;

fig = figure('Position', [10 10 900 1200])
tlo = tiledlayout(7,2,"TileSpacing","tight")
subplottcount = 1;
subplotacount = 2;

for indivIdx = 1:7
    nexttile(subplottcount)
    hold on
    nexttile(subplotacount)
    hold on

    noffsets = sum(fileIndices(:,1)==indivIdx);
    for idx = 1:length(fileIndices)
        if indivIdx
            if fileIndices(idx,1) ~= indivIdx
                continue
            end
        end
        p = struct();
        p = params_static(p);
        cutoffAngle = p.MinKneeAngle;
        p.dt = 1/2000;
        p.tfinal = 15e-4; % seconds

        p.MinKneeAngle = 0;
        
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
            force_m_k(i) = torque_m_k(i)/springArm;
            ext_m_k(i) = cos(phi(i) - p.ExtensorOffset);
            force_s_k(i) = force_m_k(i)/ext_m_k(i);
        end

        nexttile(subplottcount)
        plot(t(vq2 > cutoffAngle),-torque_m_k(vq2 > cutoffAngle),LineWidth=1.5,Color=palatte(1,:))
        plot(t(vq2 <= cutoffAngle),-torque_m_k(vq2 <= cutoffAngle),'--',LineWidth=1.5,Color=palatte(1,:))
        plot(t(vq2 > cutoffAngle),torque_m_h(vq2 > cutoffAngle),LineWidth=1.5,Color=palatte(2,:))
        plot(t(vq2 <= cutoffAngle),torque_m_h(vq2 <= cutoffAngle),'--',LineWidth=1.5,Color=palatte(2,:))
        plot(t(vq2 > cutoffAngle),torque_m_b(vq2 > cutoffAngle),LineWidth=1.5,Color=palatte(3,:))
        plot(t(vq2 <= cutoffAngle),torque_m_b(vq2 <= cutoffAngle),'--',LineWidth=1.5,Color=palatte(3,:))
        nexttile(subplotacount)
        plot(rad2deg(phi(vq2 > cutoffAngle)),-torque_m_k(vq2 > cutoffAngle),LineWidth=1.5,Color=palatte(1,:))
        plot(rad2deg(phi(vq2 <= cutoffAngle)),-torque_m_k(vq2 <= cutoffAngle),'--',LineWidth=1.5,Color=palatte(1,:))
        plot(rad2deg(phi(vq2 > cutoffAngle)),torque_m_h(vq2 > cutoffAngle),LineWidth=1.5,Color=palatte(2,:))
        plot(rad2deg(phi(vq2 <= cutoffAngle)),torque_m_h(vq2 <= cutoffAngle),'--',LineWidth=1.5,Color=palatte(2,:))
        plot(rad2deg(phi(vq2 > cutoffAngle)),torque_m_b(vq2 > cutoffAngle),LineWidth=1.5,Color=palatte(3,:))
        plot(rad2deg(phi(vq2 <= cutoffAngle)),torque_m_b(vq2 <= cutoffAngle),'--',LineWidth=1.5,Color=palatte(3,:))
    end
    if indivIdx == 7
        nexttile(subplottcount)
        ylim([-5,13]*1e-5)
        xlim([0,0.02])
        nexttile(subplotacount)
        ylim([-5,10]*1e-5)
    else
        nexttile(subplottcount)
        ylim([-5,13]*1e-4)
        xlim([0,0.045])
        nexttile(subplotacount)
        ylim([-5,10]*1e-4)
    end
    nexttile(subplottcount)
    ylabel(strcat("Individual ",num2str(indivIdx)))
    nexttile(subplotacount)
    xlim([0,110])
    xline(90,'--',LineWidth=1.5)

    subplottcount = subplottcount + 2;
    subplotacount = subplotacount + 2;
end
nexttile(2)
leg = legend(["Knee Torque" ...
    "Hip Torque" ...
    "Body Torque"],FontSize=12,Orientation="horizontal")
leg.IconColumnWidth = 10;
leg.Layout.Tile = 'north';
nexttile(13)
xlabel("$t \, \mathrm{(s)}$",Interpreter="latex",FontSize=16)
nexttile(14)
xlabel("$\theta_k \, \mathrm{(deg)}$",Interpreter="latex",FontSize=16)
ylabel(tlo,'Torque (N-m)',FontSize=16);
saveas(fig,"Figures/SpringTorqueTimeAngle.png")