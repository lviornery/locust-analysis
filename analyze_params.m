close all
clear

load("params_data.mat")

indivIdx = 3;

figure
subplot(3,2,1)
title("hhip")
hold on
subplot(3,2,2)
title("lhip")
hold on
subplot(3,2,3)
title("lbody")
hold on
subplot(3,2,4)
title("lfemur")
hold on
subplot(3,2,5)
title("ltibia")
hold on
subplot(3,2,6)
title("szr")
hold on

offset_number = 0;
meanhhip = 0;
meanlhip = 0;
meanlbody = 0;
meanlfemur = 0;
meanltibia = 0;
meanszr = 0;

noffsets = sum(fileIndices(:,1)==indivIdx);
for idx = 1:length(fileIndices)
    if indivIdx
        if fileIndices(idx,1) ~= indivIdx
            continue
        end
    end
    offset_number = offset_number + 1;
    p = struct();
    p = params_static(p);
    p.ExtensorOffset = deg2rad(5);
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
    videoName = strcat("fulldata/Locust ",num2str(fileIndices(idx,1))," jump ",num2str(fileIndices(idx,2)),".mp4");
    
    [t,vq1,vq2,vx,vy,vtheta,szr,mintthresh,maxtthresh,p] = system_extractJumpData(fileName,idx,angcutoff,p);
    
    meanhhip = meanhhip + p.hhip;
    meanlhip = meanlhip + p.lhip;
    meanlbody = meanlbody + p.lbody;
    meanlfemur = meanlfemur + p.lfemur;
    meanltibia = meanltibia + p.ltibia;
    meanszr = meanszr + szr;
    subplot(3,2,1)
    plot(offset_number,p.hhip,'o')
    subplot(3,2,2)
    plot(offset_number,p.lhip,'o')
    subplot(3,2,3)
    plot(offset_number,p.lbody,'o')
    subplot(3,2,4)
    plot(offset_number,p.lfemur,'o')
    subplot(3,2,5)
    plot(offset_number,p.ltibia,'o')
    subplot(3,2,6)
    plot(offset_number,szr,'o')
end

meanhhip = meanhhip/noffsets;
meanlhip = meanlhip/noffsets;
meanlbody = meanlbody/noffsets;
meanlfemur = meanlfemur/noffsets;
meanltibia = meanltibia/noffsets;
meanszr = meanszr/noffsets;

subplot(3,2,1)
plot([1,noffsets],[meanhhip,meanhhip],'black')
ylim(sort([meanhhip*0.6,meanhhip*1.4]))
subplot(3,2,2)
plot([1,noffsets],[meanlhip,meanlhip],'black')
ylim(sort([meanlhip*0.6,meanlhip*1.4]))
subplot(3,2,3)
plot([1,noffsets],[meanlbody,meanlbody],'black')
ylim(sort([meanlbody*0.6,meanlbody*1.4]))
subplot(3,2,4)
plot([1,noffsets],[meanlfemur,meanlfemur],'black')
ylim(sort([meanlfemur*0.6,meanlfemur*1.4]))
subplot(3,2,5)
plot([1,noffsets],[meanltibia,meanltibia],'black')
ylim(sort([meanltibia*0.6,meanltibia*1.4]))
subplot(3,2,6)
plot([1,noffsets],[meanszr,meanszr],'black')
ylim(sort([meanszr*0.6,meanszr*1.4]))

sgtitle(strcat("Individual ",num2str(indivIdx)," (",num2str(noffsets)," trajectories)"))