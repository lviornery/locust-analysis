function p = rigidbody_tform2disp(g)
%tform2disp extracts the 3-vector displacement from
%the 4x4 rigid body transform matrix g
p = g(1:3,4);
end

