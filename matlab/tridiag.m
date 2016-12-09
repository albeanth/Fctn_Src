function [ dumA, QFinal ] = tridiag( A0 )
% This function reduces any symmetric, mxm matrix to tridiagonal form
% through orthogonality similarity transforms. For symmetric cases, the
% Hessenberg form will be tridiagonal.

% It will take the original matrix, A0, and return the tridiagonal form, T,
% and associated unitary matrix Q. 

% This function should work for non-symmetric matrices too (I think). 
% It will return an upper triangular (Hessenberg) form - look into page 
% 211.

dumA = A0;
tmp = size(dumA);
m = tmp(1);
v = cell(1,length(1:m-2));

for k = 1:m-2
    k
    x = dumA(k+1:m,k);
    tmp = size(x);
    if sign(x(1)) == 0
        tmpsign = 1;
    else
        tmpsign = sign(x(1));
    end
    dum = tmpsign*norm(x,2)*eye(tmp(1),1) + x;
    if sum(x) == 0
        v{k} = dum;
        continue
    end
    v{k} = dum/norm(dum,2);
    v{k}
%     fprintf('%.15e\n',dumA(k+1:m,k:m)-2*v{k}*(v{k}'*dumA(k+1:m,k:m)))
    dumA(k+1:m,k:m) = dumA(k+1:m,k:m) - 2*v{k}*(v{k}'*dumA(k+1:m,k:m));
    2*v{k}*(v{k}'*dumA(k+1:m,k:m))
%     fprintf('l = %1.8e | %1.8e\n',dumA(1:m,k+1:m))
%     fprintf('\n')
%     fprintf('r = %1.8e | %1.8e\n',2*(dumA(1:m,k+1:m)*v{k})*v{k}')
    dumA(1:m,k+1:m) = dumA(1:m,k+1:m) - 2*(dumA(1:m,k+1:m)*v{k})*v{k}';
    2*(dumA(1:m,k+1:m)*v{k})*v{k}'
end

Q = cell(1,length(1:m-2));
QFinal = 1;
for i=1:length(v)
    l = length(v{i});
    tmp = eye(l)-2*(v{i}*v{i}');
    Itmp = eye(m);
    Itmp(i+1:m,i+1:m) = tmp;
    Q{i} = Itmp;
    QFinal = QFinal*Q{i};
end
end

