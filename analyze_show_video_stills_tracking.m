close all
clear

idx = 3;

load("params_data.mat")
load("palattes\customPalatte.mat")

p = struct();
p = params_static(p);
p.dt = 1/2000;
p.tfinal = 15e-4; % seconds

showFrame = true;

method = ["butter","movav","polyfit","tvdiff"];
method = method(4);
butterOrder = 3;
nyquistFrac = 0.2;

mmn1 = 5;
mmn2 = 9;

polyorder = 6;

a1 = 0.01;
a2 = 0.1;

cropOffset = 200
xoffset = -150
yoffset = 120

fileName = strcat("fulldata/a",num2str(fileIndices(idx,1)),"j",num2str(fileIndices(idx,2)),".csv");
videoName = strcat("fulldata/Locust ",num2str(fileIndices(idx,1))," jump ",num2str(fileIndices(idx,2)),".mp4");

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
force_m_k = zeros(length(t),1);
ext_m_k = zeros(length(t),1);
force_s_k = zeros(length(t),1);
predictForward = zeros(length(t),10);
cons_viol = zeros(length(t),4);
for i = 1:length(t)
    zi = [vq1(i);vq2(i);x(i);y(i);vtheta(i);vdq1(i);vdq2(i);dx(i);dy(i);vdtheta(i);vddq1(i);vddq2(i);ddx(i);ddy(i);vddtheta(i)];
    [torque_m_h(i),torque_m_k(i),torque_m_b(i)] = system_get_torque(zi,p);
    cons_viol(i,1:2) = gen_A(p.lfemur,p.ltibia,zi(1),zi(2),zi(5),zi(3),zi(4))*zi(6:10);
    cons_viol(i,3:4) = gen_dA(zi(6),zi(7),zi(10),zi(8),zi(9),p.lfemur,p.ltibia,zi(1),zi(2),zi(5))*zi(11:15);
    phi = pi-vq2(i);
    springArm = p.ExtensorDistance*sin(phi - p.ExtensorOffset);
    force_m_k(i) = torque_m_k(i)/springArm;
    ext_m_k(i) = cos(phi - p.ExtensorOffset);
    force_s_k(i) = force_m_k(i)/ext_m_k(i);
end

p.hipTorqueGrid = griddedInterpolant(t,torque_m_h);
p.kneeTorqueGrid = griddedInterpolant(t,torque_m_k);
p.supportTorqueGrid = griddedInterpolant(t,torque_m_b);
for i = 1:length(t)
    predictForward(i,:) = de_jump(t(i),[vq1(i),vq2(i),x(i),y(i),vtheta(i),vdq1(i),vdq2(i),dx(i),dy(i),vdtheta(i)],p);
end

%% wireframe video
p.tfinal = t(end); % seconds
p.InitialHipAngle = vq1(1);
p.InitialLegAngle = vq2(1);
p.InitialBodyAngle = vtheta(1);
p.InitialHipSpeed = vdq1(1);
p.InitialLegSpeed = vdq2(1);
p.InitialBodySpeed = vdtheta(1);

p.hipTorqueGrid = griddedInterpolant(t,torque_m_h);
p.kneeTorqueGrid = griddedInterpolant(t,torque_m_k);
p.supportTorqueGrid = griddedInterpolant(t,torque_m_b);
[tout,zmatMod,ddzout,fout,stateout] = system_calculate(p);

zmatObs = [vq1,vq2,x,y,vtheta];
zmatObs = [zmatObs, vdq1, vdq2, dx, dy, vdtheta];

for i = 1:length(t)
    positionsObs(:,:,i) = plot_jointPoints(zmatObs(i,:),p);
    positionsMod(:,:,i) = plot_jointPoints(zmatMod(i,:),p);
end

ellipsecoordinatesObs = plot_bodyellipse(zmatObs(1,:),p);
ellipsecoordinatesMod = plot_bodyellipse(zmatMod(1,:),p);

frameIndices = round(size(t,1)*linspace(0.1,0.9,5));
timeIndices = frameIndices*p.dt;
    
fig = figure('Position', [10 10 900 1200])
tlo = tiledlayout(7,5,"TileSpacing","tight")
nexttile([1,5])
colororder(palatte)
hold on
plot(t,180-rad2deg(zmatObs(:,2)),t,180-rad2deg(zmatMod(:,2)),LineWidth=1.5)
xline(timeIndices,'--',LineWidth=1.5)
hold off
ylabel("$\theta_k \, \mathrm{(deg)}$",'interpreter','latex',FontSize=14,Rotation=0)
xlim([0,t(end)])
ylim([0,100])
leg = legend(["Measured";"Feedforward"],Location="southwest",FontSize=12,Orientation="horizontal")
leg.IconColumnWidth = 10;
leg.Layout.Tile = "north"
nexttile([1,5])
colororder(palatte)
hold on
plot(t,rad2deg(zmatObs(:,1)),t,rad2deg(zmatMod(:,1)),LineWidth=1.5)
xline(timeIndices,'--',LineWidth=1.5)
hold off
ylabel("$\theta_h \, \mathrm{(deg)}$",'interpreter','latex',FontSize=14,Rotation=0)
xlim([0,t(end)])
ylim([100,180])
nexttile([1,5])
colororder(palatte)
hold on
plot(t,rad2deg(zmatObs(:,5)),t,rad2deg(zmatMod(:,5)),LineWidth=1.5)
xline(timeIndices,'--',LineWidth=1.5)
hold off
ylabel("$\theta_b \, \mathrm{(deg)}$",'interpreter','latex',FontSize=14,Rotation=0)
xlim([0,t(end)])
ylim([-25,55])
nexttile([1,5])
colororder(palatte)
hold on
plot(t,zmatObs(:,3),t,zmatMod(:,3),LineWidth=1.5)
xline(timeIndices,'--',LineWidth=1.5)
hold off
ylabel("$x \, \mathrm{(m)}$",'interpreter','latex',FontSize=14,Rotation=0)
xlim([0,t(end)])
ylim([-1e-2,3e-2])
nexttile([1,5])
colororder(palatte)
hold on
plot(t,zmatObs(:,4),t,zmatMod(:,4),LineWidth=1.5)
xline(timeIndices,'--',LineWidth=1.5)
hold off
ylabel("$y \, \mathrm{(m)}$",'interpreter','latex',FontSize=14,Rotation=0)
xlim([0,t(end)])
ylim([-1e-2,3e-2])
ax = nexttile([1,5])
colororder(palatte)
set(ax,'ColorOrderIndex',2)
hold on
plot(t,fout,LineWidth=1.5)
xline(timeIndices,'--',LineWidth=1.5)
hold off
ylabel("$\mathrm{GRF \, (N)}$",'interpreter','latex',FontSize=14,Rotation=0)
xlim([0,t(end)])

rad2deg(sqrt(mean((zmatObs(:,2)-zmatMod(:,2)).^2)))

%% embedded wreframe video
videoObject = VideoReader(videoName);
% Determine how many frames there are.
numberOfFrames = videoObject.NumFrames;
vidHeight = videoObject.Height;
vidWidth = videoObject.Width;
frameCounter = 0;
frameRange = round([(zeroOffset(2) + yoffset + [-cropOffset,cropOffset]);(-zeroOffset(1) + xoffset + [-cropOffset,cropOffset])]);
rpositions = positionsMod/szr;
rpositions(1,:,:) = -zeroOffset(1)-rpositions(1,:,:) - frameRange(2,1);
rpositions(2,:,:) = rpositions(2,:,:) + zeroOffset(2) - frameRange(1,1);
for i = 1:mintthresh
    frame = readFrame(videoObject);
end
while hasFrame(videoObject)
    frameCounter = frameCounter + 1;
    frame = readFrame(videoObject);
    if ismember(frameCounter,frameIndices)
        ax = nexttile;
        rellipsecoordinates = plot_bodyellipse(zmatMod(frameCounter,:),p);
        rellipsecoordinates = rellipsecoordinates/szr;
        rellipsecoordinates(1,:) = -zeroOffset(1) - rellipsecoordinates(1,:) - frameRange(2,1);
        rellipsecoordinates(2,:) = rellipsecoordinates(2,:) + zeroOffset(2) - frameRange(1,1);
        frame = flip(frame,1);
        frameCrop = frame(frameRange(1,1):frameRange(1,2),frameRange(2,1):frameRange(2,2),:);
        image(frameCrop,"Parent",ax);
        set(ax,'YDir','normal')
        if showFrame
            hold on
            plot([rpositions(1,4,frameCounter),rpositions(1,3,frameCounter)],[rpositions(2,4,frameCounter),rpositions(2,3,frameCounter)],LineWidth=3,Color=palatte(2,:));
            plot([rpositions(1,3,frameCounter),rpositions(1,2,frameCounter)],[rpositions(2,3,frameCounter),rpositions(2,2,frameCounter)],LineWidth=3,Color=palatte(2,:));
            plot([rpositions(1,2,frameCounter),rpositions(1,1,frameCounter)],[rpositions(2,2,frameCounter),rpositions(2,1,frameCounter)],LineWidth=3,Color=palatte(2,:));
            plot(rellipsecoordinates(1,:),rellipsecoordinates(2,:),LineWidth=1.5,Color=palatte(2,:))
            hold off
        end
        pbaspect(ax,[1 1 1])
        ax.Visible = "off";
    end
end
saveas(fig,"Figures/VideoStillsWireframe.png")