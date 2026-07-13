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
palatteskip = palatte(2:end,:);

p = struct();
p = params_static(p);
p.dt = 1/2000;
p.tfinal = 15e-4; % seconds

AllW_h = [];
AllW_k = [];
AllW_b = [];
AllW_t = [];
AllKE_o = [];
AllKE_m = [];

colorder = [2,1,5,3,4];

for indivIdx = 1:7
    for idx = 1:length(fileIndices)
        if indivIdx
            if fileIndices(idx,1) ~= indivIdx
                continue
            end
        end

        fileName = strcat("fulldata/a",num2str(fileIndices(idx,1)),"j",num2str(fileIndices(idx,2)),".csv");
        
        [t,vq1,vq2,vx,vy,vtheta,szr,mintthresh,maxtthresh,p] = system_extractJumpData(fileName,idx,p);
        
        p.MinKneeAngle = min(vq2);
        
        z0 = de_bodyICs([vq1(1);vq2(1);vtheta(1);0;0;0],p);
        
        startX = z0(3);
        startY = z0(4);
        x = (vx - vx(1))*szr+startX;
        y = (vy - vy(1))*szr+startY;
        
        zeroOffset = [vx(1) - startX/szr,vy(1) - startY/szr];
        
        %% differentiation and torque calculation
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

        torque_m_h = zeros(length(t),1);
        torque_m_k = zeros(length(t),1);
        torque_m_b = zeros(length(t),1);
        for j = 1:length(t)
            zi = [vq1(j);vq2(j);x(j);y(j);vtheta(j);vdq1(j);vdq2(j);dx(j);dy(j);vdtheta(j);vddq1(j);vddq2(j);ddx(j);ddy(j);vddtheta(j)];
            [torque_m_h(j),torque_m_k(j),torque_m_b(j)] = system_get_torque(zi,p);
        end
        
        KE_o = system_get_ke(zi,p);

        %% forward dynamics and trajectory plot
        p.tfinal = t(end); % seconds
        p.InitialHipAngle = vq1(1);
        p.InitialLegAngle = vq2(1);
        p.InitialBodyAngle = vtheta(1);
        
        p.hipTorqueGrid = griddedInterpolant(t,torque_m_h);
        p.kneeTorqueGrid = griddedInterpolant(t,torque_m_k);
        p.supportTorqueGrid = griddedInterpolant(t,torque_m_b);
        p.InitialHipSpeed = vdq1(1);
        p.InitialLegSpeed = vdq2(1);
        p.InitialBodySpeed = vdtheta(1);
        [tout,zmatMod,ddzout,fout,stateout] = system_calculate(p);

        KE_m = system_get_ke(zmatMod(end,:),p);

        P_h = torque_m_h(:).*vdq1;
        P_k = torque_m_k(:).*vdq2;
        P_b = torque_m_b(:).*vdtheta;

        W_h = cumtrapz(t,P_h);
        W_k = cumtrapz(t,P_k);
        W_b = cumtrapz(t,P_b);
        W_t = W_h + W_k + W_b;

        W_h = abs(W_h(end)/W_t(end));
        W_k = abs(W_k(end)/W_t(end));
        W_b = abs(W_b(end)/W_t(end));

        AllW_h = [AllW_h;W_h];
        AllW_k = [AllW_k;W_k];
        AllW_b = [AllW_b;W_b];
        AllW_t = [AllW_t;W_t(end)];
        AllKE_o = [AllKE_o;KE_o];
        AllKE_m = [AllKE_m;KE_m];
    end
end

fig = figure('Position',[0,0,5*150,400])
x = 1:3;
bar(x,100*[mean(AllW_k),mean(AllW_h),mean(AllW_b)],'FaceColor',palatte(2,:))
hold on
er = errorbar(x,100*[mean(AllW_k),mean(AllW_h),mean(AllW_b)],-100*[std(AllW_k),std(AllW_h),std(AllW_b)],100*[std(AllW_k),std(AllW_h),std(AllW_b)]);
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
yline(100,'--',LineWidth=2)
hold off
xticklabels({"W_k","W_h","W_b"})
ylabel("Work contributed by torque through each joint (%)")
ylim([0,105])

saveas(fig,"Figures/SI_Work.png")