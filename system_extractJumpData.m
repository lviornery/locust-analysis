function [t,vq1,vq2,vx,vy,vtheta,szr,mintthresh,maxtthresh,p] = system_extractJumpData(fileName,idx,p)

    load("params_data.mat","bodyMasses","femurMasses","tibiaMasses","ExtensorDistances","delayTimes","fileIndices");

    opts = detectImportOptions(fileName,ImportErrorRule="omitrow",MissingRule="omitrow");
    % this is import logic - depends on data format
    jump = readmatrix(fileName,opts);
    %units of time and pixels
    %format: t, comx, comy, a1, a2, a3, scalebar, hipx, hipy, flength, tlength,
    %blength
    %a1 is q2
    %a2 is theta
    %a3 is q1
    %comx and comy are the cetner of mass pixel positions
    %scalebar is the number of pixels in a 1-cm bar
    %hipx and hipy are the x and y distance from the com to the hip joint
    %flength tlength and blength (planned, currently not in use) are the pixel
    %lengths of the femur, tibia, and body
    
    % need to invert the x positions for the model to work (in model, animal
    % should jump towards the right)
    jump(:,[2,8]) = -jump(:,[2,8]);

    szr = 10*1e-3/jump(1,7); %m/px
    
    p.lfemur = jump(1,10)*szr;
    p.ltibia = jump(1,11)*szr;
    p.lbody = jump(1,12)*szr;
    
    p.mbody = bodyMasses((bodyMasses(:,1)==fileIndices(idx,1)),2);
    p.mfemur = femurMasses((femurMasses(:,1)==fileIndices(idx,1)),2);
    p.mtibia = tibiaMasses((tibiaMasses(:,1)==fileIndices(idx,1)),2);
    p.ExtensorDistance = ExtensorDistances((ExtensorDistances(:,1)==fileIndices(idx,1)),2);
    p.mtot = p.mbody + 2*p.mfemur + 2*p.mtibia;
    p.rbody = p.lbody/10; % m - body radius
    p.Ibody = 1/12*p.mbody*(3*p.rbody^2+p.lbody^2); % kg-m - moment of inertia about major axis (apumes graphopper is a cylinder)
    p.Ifemur = p.mfemur*p.lfemur^2/12;
    p.Itibia = p.mtibia*p.ltibia^2/12;
    
    vq1 = -1*jump(:,6);
    vq2 = pi + jump(:,4);
    vtheta = pi - jump(:,5);
    vx = jump(:,2);
    vy = jump(:,3);
    if vq1 > pi
        vq1 = vq1 - 2*pi;
    end
    if vtheta > pi
        vtheta = vtheta - 2*pi;
    end
    
    comhipdist = jump(1,8:9);
    comhipdist = [cos(vtheta(1)),sin(vtheta(1));-sin(vtheta(1)),cos(vtheta(1))]*comhipdist(:);
    p.lhip = comhipdist(1)*szr;
    p.hhip = comhipdist(2)*szr;
    
    mintthresh = round(delayTimes(idx)/p.dt);
    maxtthresh = find(vq2 < p.MinKneeAngle,1);
    if isempty(maxtthresh)
        maxtthresh = length(jump(:,1));
    end
    t = cumsum([0;p.dt*ones(maxtthresh-mintthresh-1,1)]);
    vq1 = vq1(mintthresh+1:maxtthresh);
    vq2 = vq2(mintthresh+1:maxtthresh);
    vx = vx(mintthresh+1:maxtthresh);
    vy = vy(mintthresh+1:maxtthresh);
    vtheta = vtheta(mintthresh+1:maxtthresh);