function w = rigidbody_skew2angvel(what)
%skew2angvel maps the 3x3 skew-symmetric
%matrix what to the 3-vector w
w = [what(3,2);what(1,3);what(2,1)];
end

