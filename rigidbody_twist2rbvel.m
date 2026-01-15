function xihat = rigidbody_twist2rbvel(xi)
%twist2rbvel maps the 6-vector xi to 
%the 4x4 matrix xihat
%initialize xihat
xihat = zeros(4);
if ~isnumeric(xi)
    xihat = sym(xihat);
end
%the upper left x3 is just what, so generate that
xihat(1:3,1:3) = rigidbody_angvel2skew(xi(4:6));
%and put v on the last column
xihat(1:3,4) = xi(1:3);
end

