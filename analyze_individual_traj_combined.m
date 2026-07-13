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

method = ["butter","movav","polyfit","tvdiff"];
method = method(4);

load("params_data.mat")
load("palattes\customPalatte.mat")

p = struct();
p = params_static(p);
p.dt = 1/2000;
p.tfinal = 15e-4; % seconds

doLineVideo = false;
doWireframe = false;
doOverlay = true;

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
end

p.hipTorqueGrid = griddedInterpolant(t,torque_m_h);
p.kneeTorqueGrid = griddedInterpolant(t,torque_m_k);
p.supportTorqueGrid = griddedInterpolant(t,torque_m_b);
for i = 1:length(t)
    predictForward(i,:) = de_jump(t(i),[vq1(i),vq2(i),x(i),y(i),vtheta(i),vdq1(i),vdq2(i),dx(i),dy(i),vdtheta(i)],p);
end
%% 
figure
subplot(1,3,1)
plot(t,vq1)
title("q1")
subplot(1,3,2)
plot(t,vdq1)
title("dq1")
subplot(1,3,3)
plot(t,vddq1)
title("ddq1")
figure
subplot(1,3,1)
plot(t,vq2)
title("q2")
subplot(1,3,2)
plot(t,vdq2)
title("dq2")
subplot(1,3,3)
plot(t,vddq2)
title("ddq2")
figure
subplot(1,3,1)
plot(t,x)
title("x")
subplot(1,3,2)
plot(t,dx)
title("dx")
subplot(1,3,3)
plot(t,ddx)
title("ddx")
figure
subplot(1,3,1)
plot(t,y)
title("y")
subplot(1,3,2)
plot(t,dy)
title("dy")
subplot(1,3,3)
plot(t,ddy)
title("ddy")
figure
subplot(1,3,1)
plot(t,vtheta)
title("theta")
subplot(1,3,2)
plot(t,vdtheta)
title("dtheta")
subplot(1,3,3)
plot(t,vddtheta)
title("ddtheta")
figure
subplot(2,3,1)
plot(t,vddq1,t,predictForward(:,6))
title("ddq1 (observed vs predicted)")
subplot(2,3,2)
plot(t,vddq2,t,predictForward(:,7))
title("ddq2 (observed vs predicted)")
subplot(2,3,3)
plot(t,ddx,t,predictForward(:,8))
title("ddx (observed vs predicted)")
subplot(2,3,4)
plot(t,ddy,t,predictForward(:,9))
title("ddy (observed vs predicted)")
subplot(2,3,5)
plot(t,vddtheta,t,predictForward(:,10))
title("ddtheta (observed vs predicted)")
figure
plot(t,cons_viol(:,1),t,cons_viol(:,2),t,cons_viol(:,3),t,cons_viol(:,4))
title("Velocity and acceleration contstraint violation")
figure
plot(t,dx,t,cumtrapz(p.dt,ddx)+dx(1))
title("Integrated acceleration and speed for each trajectory in X")
figure
plot(t,dy,t,cumtrapz(p.dt,ddy)+dy(1))
title("Integrated acceleration and speed for each trajectory in Y")
figure
plot(ext_m_k,force_m_k)
title("Spring force vs Spring Extension")
figure
subplot(1,3,1)
plot(t,torque_m_h)
title("Hip torque")
subplot(1,3,2)
plot(t,torque_m_k)
title("Knee torque")
subplot(1,3,3)
plot(t,torque_m_b)
title("Body torque")

%% forward dynamics and trajectory plot
p.tfinal = t(end); % seconds
p.InitialHipAngle = vq1(1);
p.InitialLegAngle = vq2(1);
p.InitialBodyAngle = vtheta(1);
p.InitialHipSpeed = vdq1(1);
p.InitialLegSpeed = vdq2(1);
p.InitialBodySpeed = vdtheta(1);

[tout,zmatMod,ddzout,fout,stateout] = system_calculate(p);

zmatObs = [vq1,vq2,x,y,vtheta];
% zmatObs = [zmatObs, zeros(size(zmatObs))];
zmatObs = [zmatObs, vdq1, vdq2, dx, dy, vdtheta];
titles = ["q1";"q2";"x";"y";"theta"];

figure('units','pixels','position',[0 0 2000 500])
for i = 1:5
    subplot(1,7,i)
    plot(t,zmatObs(:,i),tout,zmatMod(:,i))
    title(titles(i))
    legend(["Observed";"Modeled"])
    colororder(palatte)
end
subplot(1,7,6)
plot(tout,fout)
title("y ground reaction force")
subplot(1,7,7)
bar(tout,stateout)
title("state")
ylim([-1,4])
colororder(palatte)

%% real video
if doLineVideo
    videoObject = VideoReader(videoName);
    writerObj = VideoWriter('videos/realvideooverlay.avi');
    open(writerObj);
    % Determine how many frames there are.
    vidHeight = videoObject.Height;
    vidWidth = videoObject.Width;
    figure('Position', [10 10 vidWidth vidHeight])
    currAxes = axes;
    frameCounter = 0;
    while hasFrame(videoObject)
        frameCounter = frameCounter + 1;
        frame = readFrame(videoObject);
        if frameCounter <= size(t,1)
            frame = flip(frame ,1);
            image(frame,"Parent",currAxes);
            set(currAxes,'YDir','normal')
            hold on
            plot(1:vidWidth,(1:vidWidth)*-thighslope(frameCounter)+thighoffset(frameCounter))
            plot(1:vidWidth,(1:vidWidth)*-shinslope(frameCounter)+shinoffset(frameCounter))
            plot(currAxes,-interceptPoint(frameCounter,1),interceptPoint(frameCounter,2),'o');
            plot(1:vidWidth,(1:vidWidth)*-bodyslope(frameCounter)+bodyoffset(frameCounter))
            plot(currAxes,-x(frameCounter),y(frameCounter),'o');
            hold off
            currAxes.Visible = "off";
            F = getframe(currAxes);
            writeVideo(writerObj,F);
        end
    end
    close(writerObj);
end

%% wireframe video
for i = 1:length(t)
    positionsObs(:,:,i) = plot_jointPoints(zmatObs(i,:),p);
    positionsMod(:,:,i) = plot_jointPoints(zmatMod(i,:),p);
end
    
if doWireframe
    gcf = figure('units','pixels','position',[0 0 2000 800]);
    subplot(1,2,1)
    hold on
    ellipsecoordinatesObs = plot_bodyellipse(zmatObs(1,:),p);
    l1 = plot([positionsObs(1,4,1),positionsObs(1,3,1)],[positionsObs(2,4,1),positionsObs(2,3,1)]);
    l2 = plot([positionsObs(1,3,1),positionsObs(1,2,1)],[positionsObs(2,3,1),positionsObs(2,2,1)]);
    l3 = plot([positionsObs(1,2,1),positionsObs(1,1,1)],[positionsObs(2,2,1),positionsObs(2,1,1)]);
    l4 = plot(ellipsecoordinatesObs(1,:),ellipsecoordinatesObs(2,:));
    xlim([-.04,.04])
    ylim([-0.04,.04])
    hold off
    set(gca, 'XDir','reverse')
    title("exp")
    subplot(1,2,2)
    hold on
    ellipsecoordinatesMod = plot_bodyellipse(zmatMod(1,:),p);
    l5 = plot([positionsMod(1,4,1),positionsMod(1,3,1)],[positionsMod(2,4,1),positionsMod(2,3,1)]);
    l6 = plot([positionsMod(1,3,1),positionsMod(1,2,1)],[positionsMod(2,3,1),positionsMod(2,2,1)]);
    l7 = plot([positionsMod(1,2,1),positionsMod(1,1,1)],[positionsMod(2,2,1),positionsMod(2,1,1)]);
    l8 = plot(ellipsecoordinatesMod(1,:),ellipsecoordinatesMod(2,:));
    xlim([-.04,.04])
    ylim([-0.04,.04])
    hold off
    set(gca, 'XDir','reverse')
    title("model")
    writerObj = VideoWriter('videos/wireframevideo.avi');
    open(writerObj);
    for i = 1:length(t)
        subplot(1,2,1)
        ellipsecoordinatesObs = plot_bodyellipse(zmatObs(i,:),p);
        X1 = [positionsObs(1,4,i),positionsObs(1,3,i)];
        X2 = [positionsObs(1,3,i),positionsObs(1,2,i)];
        X3 = [positionsObs(1,2,i),positionsObs(1,1,i)];
        Y1 = [positionsObs(2,4,i),positionsObs(2,3,i)];
        Y2 = [positionsObs(2,3,i),positionsObs(2,2,i)];
        Y3 = [positionsObs(2,2,i),positionsObs(2,1,i)];
        set(l1,'XData',X1);
        set(l1,'YData',Y1);
        set(l2,'XData',X2);
        set(l2,'YData',Y2)
        set(l3,'XData',X3);
        set(l3,'YData',Y3);
        set(l4,'XData',ellipsecoordinatesObs(1,:));
        set(l4,'YData',ellipsecoordinatesObs(2,:));
    
        subplot(1,2,2)
        ellipsecoordinatesMod = plot_bodyellipse(zmatMod(i,:),p);
        X4 = [positionsMod(1,4,i),positionsMod(1,3,i)];
        X5 = [positionsMod(1,3,i),positionsMod(1,2,i)];
        X6 = [positionsMod(1,2,i),positionsMod(1,1,i)];
        Y4 = [positionsMod(2,4,i),positionsMod(2,3,i)];
        Y5 = [positionsMod(2,3,i),positionsMod(2,2,i)];
        Y6 = [positionsMod(2,2,i),positionsMod(2,1,i)];
        set(l5,'XData',X4);
        set(l5,'YData',Y4);
        set(l6,'XData',X5);
        set(l6,'YData',Y5)
        set(l7,'XData',X6);
        set(l7,'YData',Y6);
        set(l8,'XData',ellipsecoordinatesMod(1,:));
        set(l8,'YData',ellipsecoordinatesMod(2,:));
        sgtitle(['t = ', num2str(t(i),'%4.5f')]);
        F = getframe(gcf);
        writeVideo(writerObj,F);
    end
    close(writerObj);
end

%% embedded wreframe video
if doOverlay
    videoObject = VideoReader(videoName);
    writerObj = VideoWriter('videos/realvideowfoverlay.avi');
    open(writerObj);
    % Determine how many frames there are.
    numberOfFrames = videoObject.NumFrames;
    vidHeight = videoObject.Height;
    vidWidth = videoObject.Width;
    figure('Position', [10 10 vidWidth vidHeight])
    currAxes = axes;
    frameCounter = 0;
    rpositions = positionsMod/szr;
    rpositions(1,:,:) = -zeroOffset(1)-rpositions(1,:,:);
    rpositions(2,:,:) = rpositions(2,:,:) + zeroOffset(2);
    for i = 1:mintthresh
        frameCounter = 1;
        frame = readFrame(videoObject);
        frame = flip(frame,1);
        image(frame,"Parent",currAxes);
        set(currAxes,'YDir','normal')
        hold on
        rellipsecoordinates = plot_bodyellipse(zmatMod(frameCounter,:),p);
        rellipsecoordinates = rellipsecoordinates/szr;
        rellipsecoordinates(1,:) = -zeroOffset(1) - rellipsecoordinates(1,:);
        rellipsecoordinates(2,:) = rellipsecoordinates(2,:) + zeroOffset(2);
        plot([rpositions(1,4,frameCounter),rpositions(1,3,frameCounter)],[rpositions(2,4,frameCounter),rpositions(2,3,frameCounter)]);
        plot([rpositions(1,3,frameCounter),rpositions(1,2,frameCounter)],[rpositions(2,3,frameCounter),rpositions(2,2,frameCounter)]);
        plot([rpositions(1,2,frameCounter),rpositions(1,1,frameCounter)],[rpositions(2,2,frameCounter),rpositions(2,1,frameCounter)]);
        plot(rellipsecoordinates(1,:),rellipsecoordinates(2,:))
        hold off
        currAxes.Visible = "off";
        F = getframe(currAxes);
        writeVideo(writerObj,F);
    end
    framecounter = 0;
    while hasFrame(videoObject)
        frameCounter = frameCounter + 1;
        frame = readFrame(videoObject);
        if frameCounter <= size(t,1)
            frame = flip(frame,1);
            image(frame,"Parent",currAxes);
            set(currAxes,'YDir','normal')
            hold on
            rellipsecoordinates = plot_bodyellipse(zmatMod(frameCounter,:),p);
            rellipsecoordinates = rellipsecoordinates/szr;
            rellipsecoordinates(1,:) = -zeroOffset(1) - rellipsecoordinates(1,:);
            rellipsecoordinates(2,:) = rellipsecoordinates(2,:) + zeroOffset(2);
            plot([rpositions(1,4,frameCounter),rpositions(1,3,frameCounter)],[rpositions(2,4,frameCounter),rpositions(2,3,frameCounter)]);
            plot([rpositions(1,3,frameCounter),rpositions(1,2,frameCounter)],[rpositions(2,3,frameCounter),rpositions(2,2,frameCounter)]);
            plot([rpositions(1,2,frameCounter),rpositions(1,1,frameCounter)],[rpositions(2,2,frameCounter),rpositions(2,1,frameCounter)]);
            plot(rellipsecoordinates(1,:),rellipsecoordinates(2,:))
            hold off
            currAxes.Visible = "off";
            F = getframe(currAxes);
            writeVideo(writerObj,F);
        end
    end
    close(writerObj);
end