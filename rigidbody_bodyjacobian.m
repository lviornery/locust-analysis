function Jbst = rigidbody_bodyjacobian(gst,q)
%takes a symbolic homogeneous transform g and vector of generalized
%coordinates q and produces the body jacobian of the transform
Jbst = sym(zeros(6,length(q)));
invgst = inv(gst);
for i = 1:length(q)
    %the ith column is d gwb/d qi * gwb^(-1)
    Jbst(:,i) = rigidbody_rbvel2twist(invgst*diff(gst,q(i)));
end
end

