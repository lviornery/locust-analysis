function C = rigidbody_coriolismatrix(M,q,dq)
%takes a symbolic mass matrix and produces a coriolis matrix
C = sym(zeros(length(q)));
for i = 1:length(q)
    for j = 1:length(q)
        for k = 1:length(q)
            C(i,j) = C(i,j)+(diff(M(i,j),q(k))+diff(M(i,k),q(j)) - diff(M(k,j),q(i)))*dq(k)/2;
        end
    end
end
end

