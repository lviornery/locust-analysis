function Ad = rigidbody_tform2adjoint(g)
%tform2adjoint creates a 6x6 adjoint to
%the 4x4 transformation matrix g
%initialize Ad
Ad = zeros(6);
if ~isnumeric(g)
    Ad = sym(Ad);
end
%intermediate variable for the rotation matrix we pull off of g
R = g(1:3,1:3);
%R goes in the top left and bottom right
Ad(1:3,1:3) = R;
Ad(4:6,4:6) = R;
%calculate p^ (not an angular velocity, but we can use the function anyway)
%and right multiply by R
Ad(1:3,4:6) = rigidbody_angvel2skew(g(1:3,4))*R;
end

