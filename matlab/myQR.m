function [ Q,R ] = myQR( A )
% This function completes QR factorization via Householder reflectors for 
% any mxn matrix.
tmp = size(A);
m = tmp(1);
n = tmp(2);
R = A;
Q = eye(m,m);

for k=1:m-1
    x = R(k:m,k);
    if sign(x(1)) == 0
        tsign = 1;
    else
        tsign = sign(x(1));
    end
    v = tsign*norm(x,2)*eye(length(x),1) + x;
    v = v/norm(v,2);
    R(k:m,k:n) = R(k:m,k:n) - 2*v*(v'*R(k:m,k:n));
    tmpI = eye(m,m);
    tmpI(k:m,k:m) = tmpI(k:m,k:m)-2*(v*v');
    Q = Q*tmpI;
end