function [ soln ] = QRsolve( A,b )
% This function uses QR factoriation of A to solve a system of equations Ax=b.
% the solution to Ax=b, x, is returned.
tmp = size(A);
m = tmp(1);
n = tmp(2);
R = A;

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
    b(k:m) = b(k:m) - 2*v*(v'*b(k:m));    
end
soln = backsub(R,b);

function x = backsub(U,y)
% Solve U*x=y by backsubstitution
% U must be upper triangular and have at least as many rows as columns. 

[m,n]=size(U);
if (length(y) ~= m) | (m < n) 
    disp(['Error in backsub.m when solving Ux=y!'])
    disp(['Invalid matrix sizes of U and/or y'])
    sizeU = size(U)
    sizey = size(y)
    return
end

du = diag(U);

if nnz(du) < n
    disp(['Error in backsub.m !'])
    disp([' Matrix U has a zero on the diagonal.'])
    return
end

x = zeros(n,1);
yc = zeros(m,1); yc(1:m) = y(1:m);
for k = n:-1:1
    x(k) = (yc(k) - U(k,k+1:n)*x(k+1:n))/du(k);
end

