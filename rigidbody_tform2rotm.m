function R = rigidbody_tform2rotm(g)
%tform2rotm extracts the 3x3 rotation matrix from
%the 4x4 rigid body transform matrix g
R = g(1:3,1:3);
end