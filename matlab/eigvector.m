function [ v ] = eigvector( eigvals,A0,m )
% This function uses inverse iteration to obtain the eigenvectors of known
% eigenvalues.

warning('off','all')
for val=eigvals
    dumv = eye(m,1);
    for k=1:3
        w = (A0-val*eye(m,m))\dumv;
        dumv = w/norm(w,2);
    end
    if val == eigvals(1)
        v = [dumv];
    else
        v = [v dumv];
    end

end

end