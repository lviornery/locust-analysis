close all
clear

butterOrder = 3;
nyquistFrac = 0.2;

mmn1 = 5;
mmn2 = 9;

polyorder = 6;

a1 = 0.1;
a2 = 0.1;

load("params_data.mat")
load("palattes\customPalatte.mat")
palatteskip = palatte(2:end,:);

p = struct();
p = params_static(p);
p.dt = 1/2000;
p.tfinal = 15e-4; % seconds

AllRMS_b = [];
AllRMS_m = [];
AllRMS_p = [];
AllRMS_t = [];

colorder = [2,1,5,3,4];

violinorder = [1,3,4,5];

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

        %butter
        [b,a] = butter(butterOrder,nyquistFrac);
        vq1_b = filtfilt(b,a,vq1);
        vq2_b = filtfilt(b,a,vq2);
        x_b = filtfilt(b,a,x);
        y_b = filtfilt(b,a,y);
        vtheta_b = filtfilt(b,a,vtheta);
    
        vdq1_b = Diff*vq1_b;
        vdq2_b = Diff*vq2_b;
        dx_b = Diff*x_b;
        dy_b = Diff*y_b;
        vdtheta_b = Diff*vtheta_b;
        
        vddq1_b = Diff2*vq1_b;
        vddq2_b = Diff2*vq2_b;
        ddx_b = Diff2*x_b;
        ddy_b = Diff2*y_b;
        vddtheta_b = Diff2*vtheta_b;

        torque_m_h_b = zeros(length(t),1);
        torque_m_k_b = zeros(length(t),1);
        torque_m_b_b = zeros(length(t),1);
        cons_viol_b = zeros(length(t),4);
        for j = 1:length(t)
            zi = [vq1(j);vq2(j);x(j);y(j);vtheta(j);vdq1_b(j);vdq2_b(j);dx_b(j);dy_b(j);vdtheta_b(j);vddq1_b(j);vddq2_b(j);ddx_b(j);ddy_b(j);vddtheta_b(j)];
            [torque_m_h_b(j),torque_m_k_b(j),torque_m_b_b(j)] = system_get_torque(zi,p);
            cons_viol_b(j,1:2) = gen_A(p.lfemur,p.ltibia,zi(1),zi(2),zi(5),zi(3),zi(4))*zi(6:10);
            cons_viol_b(j,3:4) = gen_dA(zi(6),zi(7),zi(10),zi(8),zi(9),p.lfemur,p.ltibia,zi(1),zi(2),zi(5))*zi(11:15);
        end
        cons_viol_b = max(abs(cons_viol_b'));
        
        %movav  
        vdq1_m = movmean(Diff*vq1,mmn1);
        vdq2_m = movmean(Diff*vq2,mmn1);
        dx_m = movmean(Diff*x,mmn1);
        dy_m = movmean(Diff*y,mmn1);
        vdtheta_m = movmean(Diff*vtheta,mmn1);
        
        vddq1_m = movmean(Diff2*vq1,mmn2);
        vddq2_m = movmean(Diff2*vq2,mmn2);
        ddx_m = movmean(Diff2*x,mmn2);
        ddy_m = movmean(Diff2*y,mmn2);
        vddtheta_m = movmean(Diff2*vtheta,mmn2);

        torque_m_h_m = zeros(length(t),1);
        torque_m_k_m = zeros(length(t),1);
        torque_m_b_m = zeros(length(t),1);
        cons_viol_m = zeros(length(t),4);
        for j = 1:length(t)
            zi = [vq1(j);vq2(j);x(j);y(j);vtheta(j);vdq1_m(j);vdq2_m(j);dx_m(j);dy_m(j);vdtheta_m(j);vddq1_m(j);vddq2_m(j);ddx_m(j);ddy_m(j);vddtheta_m(j)];
            [torque_m_h_m(j),torque_m_k_m(j),torque_m_b_m(j)] = system_get_torque(zi,p);
            cons_viol_m(j,1:2) = gen_A(p.lfemur,p.ltibia,zi(1),zi(2),zi(5),zi(3),zi(4))*zi(6:10);
            cons_viol_m(j,3:4) = gen_dA(zi(6),zi(7),zi(10),zi(8),zi(9),p.lfemur,p.ltibia,zi(1),zi(2),zi(5))*zi(11:15);
        end
        cons_viol_m = max(abs(cons_viol_m'));
        
        %polyfit
        pq1 = polyfit(t,vq1,polyorder);
        pq2 = polyfit(t,vq2,polyorder);
        px = polyfit(t,x,polyorder);
        py = polyfit(t,y,polyorder);
        ptheta = polyfit(t,vtheta,polyorder);
        
        vdq1_p = polyval(polyder(pq1),t);
        vdq2_p = polyval(polyder(pq2),t);
        dx_p = polyval(polyder(px),t);
        dy_p = polyval(polyder(py),t);
        vdtheta_p = polyval(polyder(ptheta),t);
        
        vddq1_p = polyval(polyder(polyder(pq1)),t);
        vddq2_p = polyval(polyder(polyder(pq2)),t);
        ddx_p = polyval(polyder(polyder(px)),t);
        ddy_p = polyval(polyder(polyder(py)),t);
        vddtheta_p = polyval(polyder(polyder(ptheta)),t);

        torque_m_h_p = zeros(length(t),1);
        torque_m_k_p = zeros(length(t),1);
        torque_m_b_p = zeros(length(t),1);
        cons_viol_p = zeros(length(t),4);
        for j = 1:length(t)
            zi = [vq1(j);vq2(j);x(j);y(j);vtheta(j);vdq1_p(j);vdq2_p(j);dx_p(j);dy_p(j);vdtheta_p(j);vddq1_p(j);vddq2_p(j);ddx_p(j);ddy_p(j);vddtheta_p(j)];
            [torque_m_h_p(j),torque_m_k_p(j),torque_m_b_p(j)] = system_get_torque(zi,p);
            cons_viol_p(j,1:2) = gen_A(p.lfemur,p.ltibia,zi(1),zi(2),zi(5),zi(3),zi(4))*zi(6:10);
            cons_viol_p(j,3:4) = gen_dA(zi(6),zi(7),zi(10),zi(8),zi(9),p.lfemur,p.ltibia,zi(1),zi(2),zi(5))*zi(11:15);
        end
        cons_viol_p = max(abs(cons_viol_p'));

        %tvdiff
        [vdq1_t,vddq1_t] = regDiff(vq1,a1,a2,p.dt,[],1e10,1e-6,1e-10,5e-8);
        [vdq2_t,vddq2_t] = regDiff(vq2,a1,a2,p.dt,[],1e10,1e-6,1e-10,5e-8);
        [dx_t,ddx_t] = regDiff(x,a1,a2,p.dt,[],1e10,1e-6,1e-10,5e-8);
        [dy_t,ddy_t] = regDiff(y,a1,a2,p.dt,[],1e10,1e-6,1e-10,5e-8);
        [vdtheta_t,vddtheta_t] = regDiff(vtheta,a1,a2,p.dt,[],1e10,1e-6,1e-10,5e-8);

        torque_m_h_t = zeros(length(t),1);
        torque_m_k_t = zeros(length(t),1);
        torque_m_b_t = zeros(length(t),1);
        cons_viol_t = zeros(length(t),4);
        for j = 1:length(t)
            zi = [vq1(j);vq2(j);x(j);y(j);vtheta(j);vdq1_t(j);vdq2_t(j);dx_t(j);dy_t(j);vdtheta_t(j);vddq1_t(j);vddq2_t(j);ddx_t(j);ddy_t(j);vddtheta_t(j)];
            [torque_m_h_t(j),torque_m_k_t(j),torque_m_b_t(j)] = system_get_torque(zi,p);
            cons_viol_t(j,1:2) = gen_A(p.lfemur,p.ltibia,zi(1),zi(2),zi(5),zi(3),zi(4))*zi(6:10);
            cons_viol_t(j,3:4) = gen_dA(zi(6),zi(7),zi(10),zi(8),zi(9),p.lfemur,p.ltibia,zi(1),zi(2),zi(5))*zi(11:15);
        end
        cons_viol_t = max(abs(cons_viol_t'));

        %% forward dynamics and trajectory plot
        p.tfinal = t(end); % seconds
        p.InitialHipAngle = vq1(1);
        p.InitialLegAngle = vq2(1);
        p.InitialBodyAngle = vtheta(1);

        %butter
        p.hipTorqueGrid = griddedInterpolant(t,torque_m_h_b);
        p.kneeTorqueGrid = griddedInterpolant(t,torque_m_k_b);
        p.supportTorqueGrid = griddedInterpolant(t,torque_m_b_b);
        p.InitialHipSpeed = vdq1_b(1);
        p.InitialLegSpeed = vdq2_b(1);
        p.InitialBodySpeed = vdtheta_b(1);
        [tout_b,zmatMod_b,ddzout_b,fout_b,stateout_b] = system_calculate(p);
        transitionIndex_b = find(stateout_b ~= stateout_b(1),1);

        %movav
        p.hipTorqueGrid = griddedInterpolant(t,torque_m_h_m);
        p.kneeTorqueGrid = griddedInterpolant(t,torque_m_k_m);
        p.supportTorqueGrid = griddedInterpolant(t,torque_m_b_m);
        p.InitialHipSpeed = vdq1_m(1);
        p.InitialLegSpeed = vdq2_m(1);
        p.InitialBodySpeed = vdtheta_m(1);
        [tout_m,zmatMod_m,ddzout_m,fout_m,stateout_m] = system_calculate(p);
        transitionIndex_m = find(stateout_m ~= stateout_m(1),1);

        %polyfit
        p.hipTorqueGrid = griddedInterpolant(t,torque_m_h_p);
        p.kneeTorqueGrid = griddedInterpolant(t,torque_m_k_p);
        p.supportTorqueGrid = griddedInterpolant(t,torque_m_b_p);
        p.InitialHipSpeed = vdq1_p(1);
        p.InitialLegSpeed = vdq2_p(1);
        p.InitialBodySpeed = vdtheta_p(1);
        [tout_p,zmatMod_p,ddzout_p,fout_p,stateout_p] = system_calculate(p);
        transitionIndex_p = find(stateout_p ~= stateout_p(1),1);

        %tvdiff
        p.hipTorqueGrid = griddedInterpolant(t,torque_m_h_t);
        p.kneeTorqueGrid = griddedInterpolant(t,torque_m_k_t);
        p.supportTorqueGrid = griddedInterpolant(t,torque_m_b_t);
        p.InitialHipSpeed = vdq1_t(1);
        p.InitialLegSpeed = vdq2_t(1);
        p.InitialBodySpeed = vdtheta_t(1);
        [tout_t,zmatMod_t,ddzout_t,fout_t,stateout_t] = system_calculate(p);
        transitionIndex_t = find(stateout_t ~= stateout_t(1),1);
        
        zmatObs = [vq1,vq2,x,y,vtheta];

        %RMSEs
        RMS_b = rmse(zmatObs(:,colorder),zmatMod_b(:,colorder));
        RMS_m = rmse(zmatObs(:,colorder),zmatMod_m(:,colorder));
        RMS_p = rmse(zmatObs(:,colorder),zmatMod_p(:,colorder));
        RMS_t = rmse(zmatObs(:,colorder),zmatMod_t(:,colorder));
        AllRMS_b = [AllRMS_b;RMS_b];
        AllRMS_m = [AllRMS_m;RMS_m];
        AllRMS_p = [AllRMS_p;RMS_p];
        AllRMS_t = [AllRMS_t;RMS_t];
        RMS_names = ["$\theta_k RMS$";"$\theta_h RMS$";"$\theta_b RMS$";"$x RMS$";"$y RMS$"];

        RMStable = table(RMS_names,RMS_b',RMS_m',RMS_p',RMS_t',...
            'VariableNames',["Var";"B";"MA";"P";"TV"],...
            'RowNames',["q1RMS","q2RMS","xRMS","yRMS","thetaRMS",]);
        titles = {'$\theta_h$ (deg)';'$\theta_k$ (deg)';'$x$ (m)';'$y$ (m)';'$\theta_b$ (deg)';'State';...
            '$\dot{\theta_h}$ (deg/s)';'$\dot{\theta_k}$ (deg/s)';'$\dot{x}$ (m/s)';'$\dot{y}$ (m/s)';'$\dot{\theta_b}$ (deg/s)';'Constraint Violation';...
            '$\ddot{\theta_h}$ (deg/s$^2$)';'$\ddot{\theta_k}$ (deg/s$^2$)';'$\ddot{x}$ (m/s$^2$)';'$\ddot{y}$ (m/s$^2$)';'$\ddot{\theta_b}$ (deg/s$^2$)';'Ground reaction force'};
        
        fig = figure('units','pixels','position',[0 0 2000 3500],'color','white');
        axarray = [];
        for i = 1:5
            j = colorder(i);
            if j == 1 || j == 2 || j == 3
                thisZMatObs = rad2deg(zmatObs(:,j));
                thisZMatModB = rad2deg([zmatMod_b(:,j),zmatMod_b(:,j+5),ddzout_b(:,j)]);
                thisZMatModM = rad2deg([zmatMod_m(:,j),zmatMod_m(:,j+5),ddzout_m(:,j)]);
                thisZMatModP = rad2deg([zmatMod_p(:,j),zmatMod_p(:,j+5),ddzout_p(:,j)]);
                thisZMatModT = rad2deg([zmatMod_t(:,j),zmatMod_t(:,j+5),ddzout_t(:,j)]);
            end
            axarray = [axarray,subplot(3,6,i)];
            plot(t,thisZMatObs,tout_b,thisZMatModB(:,1),tout_m,thisZMatModM(:,1),tout_p,thisZMatModP(:,1),tout_t,thisZMatModT(:,1))
            ylabel(titles(j),'Interpreter','latex',FontSize=16)
            legend(["Observed";"Butterworth";"Moving Average";"Polynomial";"TV-Diff"],Location="northwest")
            subplot(3,6,i+6)
            plot(tout_b,thisZMatModB(:,2),tout_m,thisZMatModM(:,2),tout_p,thisZMatModP(:,2),tout_t,thisZMatModT(:,2))
            ylabel(titles(j+6),'Interpreter','latex',FontSize=16)
            legend(["Butterworth";"Moving Average";"Polynomial";"TV-Diff"])
            subplot(3,6,i+12)
            plot(tout_b,thisZMatModB(:,3),tout_m,thisZMatModM(:,3),tout_p,thisZMatModP(:,3),tout_t,thisZMatModT(:,3))
            ylabel(titles(j+12),'Interpreter','latex',FontSize=16)
            xlabel("t (s)")
            legend(["Butterworth";"Moving Average";"Polynomial";"TV-Diff"])
        end
        ha = subplot(3,6,6);
        pos = get(ha,'Position');
        un = get(ha,'Units');
        delete(ha)
        ht = uitable("RowName",[],"Data",RMStable,'ColumnWidth',{40.5,40.5,40.5,40.5,40.5},'Units',un,'Position',pos);
        s = uistyle("Interpreter","latex");
        addStyle(ht,s)

        subplot(3,6,12)
        plot(t,cons_viol_b,t,cons_viol_m,t,cons_viol_p,t,cons_viol_t)
        legend(["Butterworth";"Moving Average";"Polynomial";"TV-Diff"])
        ylabel(titles(12),'Interpreter','latex',FontSize=16)
        normalXLim = xlim;

        subplot(3,6,18)
        plot(tout_b,fout_b,tout_m,fout_m,tout_p,fout_p,tout_t,fout_t)
        legend(["Butterworth";"Moving Average";"Polynomial";"TV-Diff"])
        ylabel(titles(18),'Interpreter','latex',FontSize=16)
        xlabel("t (s)")

        colororder(palatteskip)
        for j = 1:5
            colororder(axarray(j),palatte)
        end
        sgtitle(strcat("Animal ",num2str(fileIndices(idx,1)),", jump ",num2str(fileIndices(idx,2))))

        exportapp(fig,strcat("Figures/indivTraj/a",num2str(fileIndices(idx,1)),"j",num2str(fileIndices(idx,2)),".png"))
        close all
    end
end



%% tbhpropaganda, buillt II thinkre
fig = figure('Position',[0,0,5*150,400]);
titles = {"$RMS_{\theta_k}$ (deg)";"$RMS_{\theta_h}$ (deg)";"$RMS_{\theta_b}$ (deg)";"$RMS_{x}$ (m)";"$RMS_{y}$ (m)"};
for j = 1:5
    AllRMS_var = [AllRMS_b(:,j),AllRMS_m(:,j),AllRMS_p(:,j),AllRMS_t(:,j)];
    if j == 1 || j == 2 || j == 3
        AllRMS_var = rad2deg(AllRMS_var);
    end
    AvgRMS_var = mean(AllRMS_var);
    ax = subplot(1,5,j);
    temp_line = plot(categorical(1:4),zeros(1,4));
    delete(temp_line)
    hold on
    for k = 1:4
        if j == 4|| j == 5
            [f,xf] = kde(AllRMS_var(:,k),Kernel="normal",Bandwidth=0.00015,Support="nonnegative",NumPoints=1000);
        else
            [f,xf] = kde(AllRMS_var(:,k),Kernel="normal",Bandwidth=0.5,Support="nonnegative",NumPoints=1000);
        end
        vplot = violinplot(EvaluationPoints=xf,DensityValues=f,FaceColorMode="manual",FaceColor=palatte(violinorder(k),:));
        vplot.XData = k;
        vplot.SeriesIndex = k;
        density_average = interp1(vplot.EvaluationPoints,...
                             vplot.DensityValues,...
                             AvgRMS_var(k));
        width_average = vplot.DensityWidth.*density_average./max(vplot.DensityValues);
        x_values = k + [-0.5 0.5]*width_average;
        y_values = [AvgRMS_var(k) AvgRMS_var(k)];
        plot(x_values, y_values, SeriesIndex=vplot.SeriesIndex,LineWidth=2,Color=palatte(violinorder(k),:));
        if j == 4|| j == 5
            ylim([0,0.01]);
        else
            ylim([0,45]);
        end
    end
    hold off
    xticklabels({"Butterworth";"Moving Average";"Polynomial";"TV-Diff"})
    title(titles(j),'Interpreter','latex',FontSize=10)
    colororder(palatteskip)

    saveas(fig,"Figures/Figure_5.png")
end