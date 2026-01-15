function xi = rigidbody_rbvel2twist(xihat)
%twist2rbvel maps the 4x4 matrix xihat to
%the 6-vector xi
xi = [xihat(1:3,4);rigidbody_skew2angvel(xihat(1:3,1:3))];
end

