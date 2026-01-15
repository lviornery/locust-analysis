function Jsst = rigidbody_spatialjacobian(gst,q)
%takes a symbolic homogeneous transform g and vector of generalized
%coordinates q and produces the spatial jacobian of the transform
Jsst = sym(zeros(6,length(q)));
invgst = inv(gst);
for i = 1:length(q)
    %the ith column is d gwb/d qi * gwb^(-1)
    Jsst(:,i) = rigidbody_rbvel2twist(diff(gst,q(i))*invgst);
end
end