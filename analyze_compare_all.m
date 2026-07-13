close all
clear

butterOrder = 3;
nyquistFrac = 0.2;

mmn1 = 5;
mmn2 = 9;

polyorder = 6;

a1 = 0.1;
a2 = 0.1;

idx = 3;

load("params_data.mat")
load("palattes\customPalatte.mat")
palatteskip = palatte(2:end,:);

p = struct();
p = params_static(p);
p.dt = 1/2000;
p.tfinal = 15e-4; % seconds

fileName = strcat("fulldata/a",num2str(fileIndices(idx,1)),"j",num2str(fileIndices(idx,2)),".csv");

[t,vq1,vq2,vx,vy,vtheta,szr,mintthresh,maxtthresh,p] = system_extractJumpData(fileName,idx,p);

z0 = de_bodyICs([vq1(1);vq2(1);vtheta(1);0;0;0],p);

startX = z0(3);
startY = z0(4);
x = (vx - vx(1))*szr+startX;
y = (vy - vy(1))*szr+startY;

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
for i = 1:length(t)
    zi = [vq1(i);vq2(i);x(i);y(i);vtheta(i);vdq1_b(i);vdq2_b(i);dx_b(i);dy_b(i);vdtheta_b(i);vddq1_b(i);vddq2_b(i);ddx_b(i);ddy_b(i);vddtheta_b(i)];
    [torque_m_h_b(i),torque_m_k_b(i),torque_m_b_b(i)] = system_get_torque(zi,p);
    cons_viol_b(i,1:2) = gen_A(p.lfemur,p.ltibia,zi(1),zi(2),zi(5),zi(3),zi(4))*zi(6:10);
    cons_viol_b(i,3:4) = gen_dA(zi(6),zi(7),zi(10),zi(8),zi(9),p.lfemur,p.ltibia,zi(1),zi(2),zi(5))*zi(11:15);
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
for i = 1:length(t)
    zi = [vq1(i);vq2(i);x(i);y(i);vtheta(i);vdq1_m(i);vdq2_m(i);dx_m(i);dy_m(i);vdtheta_m(i);vddq1_m(i);vddq2_m(i);ddx_m(i);ddy_m(i);vddtheta_m(i)];
    [torque_m_h_m(i),torque_m_k_m(i),torque_m_b_m(i)] = system_get_torque(zi,p);
    cons_viol_m(i,1:2) = gen_A(p.lfemur,p.ltibia,zi(1),zi(2),zi(5),zi(3),zi(4))*zi(6:10);
    cons_viol_m(i,3:4) = gen_dA(zi(6),zi(7),zi(10),zi(8),zi(9),p.lfemur,p.ltibia,zi(1),zi(2),zi(5))*zi(11:15);
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
for i = 1:length(t)
    zi = [vq1(i);vq2(i);x(i);y(i);vtheta(i);vdq1_p(i);vdq2_p(i);dx_p(i);dy_p(i);vdtheta_p(i);vddq1_p(i);vddq2_p(i);ddx_p(i);ddy_p(i);vddtheta_p(i)];
    [torque_m_h_p(i),torque_m_k_p(i),torque_m_b_p(i)] = system_get_torque(zi,p);
    cons_viol_p(i,1:2) = gen_A(p.lfemur,p.ltibia,zi(1),zi(2),zi(5),zi(3),zi(4))*zi(6:10);
    cons_viol_p(i,3:4) = gen_dA(zi(6),zi(7),zi(10),zi(8),zi(9),p.lfemur,p.ltibia,zi(1),zi(2),zi(5))*zi(11:15);
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
for i = 1:length(t)
    zi = [vq1(i);vq2(i);x(i);y(i);vtheta(i);vdq1_t(i);vdq2_t(i);dx_t(i);dy_t(i);vdtheta_t(i);vddq1_t(i);vddq2_t(i);ddx_t(i);ddy_t(i);vddtheta_t(i)];
    [torque_m_h_t(i),torque_m_k_t(i),torque_m_b_t(i)] = system_get_torque(zi,p);
    cons_viol_t(i,1:2) = gen_A(p.lfemur,p.ltibia,zi(1),zi(2),zi(5),zi(3),zi(4))*zi(6:10);
    cons_viol_t(i,3:4) = gen_dA(zi(6),zi(7),zi(10),zi(8),zi(9),p.lfemur,p.ltibia,zi(1),zi(2),zi(5))*zi(11:15);
end
cons_viol_t = max(abs(cons_viol_t'));

fig = figure(Position=[10,10,5*150,600])
tlo = tiledlayout(5,1,"TileSpacing","tight")
tlk = tiledlayout(tlo,1,3,"TileSpacing","tight")
tlk.Layout.Tile = 1;
tlh = tiledlayout(tlo,1,3,"TileSpacing","tight")
tlh.Layout.Tile = 2;
tlb = tiledlayout(tlo,1,3,"TileSpacing","tight")
tlb.Layout.Tile = 3;
tlx = tiledlayout(tlo,1,3,"TileSpacing","tight")
tlx.Layout.Tile = 4;
tly = tiledlayout(tlo,1,3,"TileSpacing","tight")
tly.Layout.Tile = 5;

nexttile(tlk)
plot(t,rad2deg(pi - vq2),LineWidth=1.5)
ylabel("$\theta_k \, (\textrm{deg})$",interpreter="latex",FontSize=10)
title("Position",'interpreter','latex',FontSize=10)
colororder(palatte)
nexttile(tlk)
plot([NaN],[NaN],t,rad2deg(-vdq2_m),t,rad2deg(-vdq2_p),t,rad2deg(-vdq2_b),t,rad2deg(-vdq2_t),LineWidth=1.5)
ylabel("$\dot{\theta}_k \, (\textrm{deg/s})$",interpreter="latex",FontSize=10)
title("Velocity",'interpreter','latex',FontSize=10)
colororder(palatte)
nexttile(tlk)
plot([NaN],[NaN],t,rad2deg(-vddq2_m),t,rad2deg(-vddq2_p),t,rad2deg(-vddq2_b),t,rad2deg(-vddq2_t),LineWidth=1.5)
ylabel("$\ddot{\theta}_k \, (\textrm{deg/s}^2)$",interpreter="latex",FontSize=10)
title("Acceleration",'interpreter','latex',FontSize=10)
ylim([-0.5,5.5]*1e5)
colororder(palatte)

nexttile(tlh)
plot(t,rad2deg(vq1),LineWidth=1.5)
ylabel("$\theta_h \, (\textrm{deg})$",interpreter="latex",FontSize=10)
colororder(palatte)
nexttile(tlh)
plot([NaN],[NaN],t,rad2deg(vdq1_m),t,rad2deg(vdq1_p),t,rad2deg(vdq1_b),t,rad2deg(vdq1_t),LineWidth=1.5)
ylabel("$\dot{\theta}_h \, (\textrm{deg/s})$",interpreter="latex",FontSize=10)
colororder(palatte)
nexttile(tlh)
plot([NaN],[NaN],t,rad2deg(vddq1_m),t,rad2deg(vddq1_p),t,rad2deg(vddq1_b),t,rad2deg(vddq1_t),LineWidth=1.5)
ylabel("$\ddot{\theta}_h \, (\textrm{deg/s}^2)$",interpreter="latex",FontSize=10)
ylim([-1,4]*1e5)
colororder(palatte)

nexttile(tlb)
plot(t,rad2deg(vtheta),LineWidth=1.5)
ylabel("$\theta_b \, (\textrm{deg})$",interpreter="latex",FontSize=10)
colororder(palatte)
nexttile(tlb)
plot([NaN],[NaN],t,rad2deg(vdtheta_m),t,rad2deg(vdtheta_p),t,rad2deg(vdtheta_b),t,rad2deg(vdtheta_t),LineWidth=1.5)
ylabel("$\dot{\theta}_b \, (\textrm{deg/s})$",interpreter="latex",FontSize=10)
colororder(palatte)
nexttile(tlb)
plot([NaN],[NaN],t,rad2deg(vddtheta_m),t,rad2deg(vddtheta_p),t,rad2deg(vddtheta_b),t,rad2deg(vddtheta_t),LineWidth=1.5)
ylabel("$\ddot{\theta}_b \, (\textrm{deg/s}^2)$",interpreter="latex",FontSize=10)
ylim([-1.5,1.5]*1e5)
colororder(palatte)

nexttile(tlx)
plot(t,x,LineWidth=1.5)
ylabel("$x \, (\textrm{m})$",interpreter="latex",FontSize=10)
nexttile(tlx)
plot([NaN],[NaN],t,dx_m,t,dx_p,t,dx_b,t,dx_t,LineWidth=1.5)
ylabel("$\dot{x} \, (\textrm{m/s})$",interpreter="latex",FontSize=10)
nexttile(tlx)
plot([NaN],[NaN],t,ddx_m,t,ddx_p,t,ddx_b,t,ddx_t,LineWidth=1.5)
ylabel("$\ddot{x} \, (\textrm{m/s}^2)$",interpreter="latex",FontSize=10)
colororder(palatte)

nexttile(tly)
plot(t,y,LineWidth=1.5)
ylabel("$y \, (\textrm{m})$",interpreter="latex",FontSize=10)
nexttile(tly)
plot([NaN],[NaN],t,dy_m,t,dy_p,t,dy_b,t,dy_t,LineWidth=1.5)
ylabel("$\dot{y} \, (\textrm{m/s})$",interpreter="latex",FontSize=10)
% ylim(1.25*[min(dy-eq_dy),max(dy+eq_dy)])
nexttile(tly)
plot([NaN],[NaN],t,ddy_m,t,ddy_p,t,ddy_b,t,ddy_t,LineWidth=1.5)
ylabel("$\ddot{y} \, (\textrm{m/s}^2)$",interpreter="latex",FontSize=10)
colororder(palatte)

xlabel(tlo,'Time (s)',FontSize=10);
ax = nexttile(tlk,2)
leg = legend(...
    [ax.Children(5),ax.Children(4),ax.Children(3),ax.Children(2),ax.Children(1)],...
    ["Measured";"Moving-Average Filtered";"Polynomial Fit";"Butterworth Filtered";"TV Differentiated"],...
    FontSize=10,Orientation="horizontal")
leg.IconColumnWidth = 10;
leg.Layout.Tile = 'north';
annotation("textbox",String="A",FontSize=16,LineStyle="none",Position=[0.025,0.9,0.1,0.1],FontName='Palatino Linotype')

saveas(fig,"Figures/Figure_3_A.png")

fig = figure(Position=[10,10,5*150,200])
tlo = tiledlayout(1,3,"TileSpacing","tight")
nexttile
plot([NaN],[NaN],t,-torque_m_k_m,t,-torque_m_k_p,t,-torque_m_k_b,t,-torque_m_k_t,LineWidth=1.5)
title("$\tau_k$",'interpreter','latex',FontSize=10)
ylim([-1,7]*1e-4)
% ylim(1.25*[min(torque_m_k-eq_torque_k),max(torque_m_k+eq_torque_k)])
nexttile
plot([NaN],[NaN],t,torque_m_h_m,t,torque_m_h_p,t,torque_m_h_b,t,torque_m_h_t,LineWidth=1.5)
title("$\tau_h$",'interpreter','latex',FontSize=10)
ylim([-1.5,1.5]*1e-4)
% ylim(1.25*[min(torque_m_h-eq_torque_h),max(torque_m_h+eq_torque_h)])
ax = nexttile;
plot([NaN],[NaN],t,torque_m_b_m,t,torque_m_b_p,t,torque_m_b_b,t,torque_m_b_t,LineWidth=1.5)
title("$\tau_b$",'interpreter','latex',FontSize=10)
ylim([-1.5,2.5]*1e-4)
% ylim(1.25*[min(torque_m_b-eq_torque_b),max(torque_m_b+eq_torque_b)])
colororder(palatte)
ylabel(tlo,'Torque (N-m)',FontSize=10);
xlabel(tlo,'Time (s)',FontSize=10);
leg = legend(...
    [ax.Children(4),ax.Children(3),ax.Children(2),ax.Children(1)],...
    ["Moving-Average Filtered";"Polynomial Fit";"Butterworth Filtered";"TV Differentiated"],...
    FontSize=10,Orientation="horizontal")
leg.IconColumnWidth = 10;
leg.Layout.Tile = 'north';
annotation("textbox",String="B",FontSize=16,LineStyle="none",Position=[0.025,0.9,0.1,0.1],FontName='Palatino Linotype')

saveas(fig,"Figures/Figure_3_B.png")