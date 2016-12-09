clear,clc, close
format long
A0 = [1,2,3,4;5,6,7,8;9,10,11,12;13,14,15,16];
% A0 = hilb(4);
% A0 = [0,1,0,0;1,0,0,0;0,0,0,1;0,0,1,0];
% A0 = diag(15:-1:1) + ones(15,15);
m = length(A0);

%% Part C - (Driver part)
[ H, Q ] = tridiag( A0 );
flag = 1; % running Wilkinson shifted QR algorith.
% flag = 0; % running unshifted QR algorith.
EigInfo = cell(1,m);
% fprintf('Eigenvalues of A0 are:\n')
plt = [];
eigv = [];
for i = 1:m-1
    [ H, epsi ] = qralg( H(1:m+1-i,1:m+1-i), flag );
    EigInfo{i} = [H(end,end), epsi];
%     fprintf('%.12e\n',EigInfo{i}(1))
    eigv = [eigv EigInfo{i}(1)];
    plt = [plt EigInfo{i}(2:end)];
    if length(H) == 2
        EigInfo{i+1} = [H(end-1,end-1), epsi];
%         fprintf('%.12e\n',EigInfo{i+1}(1))
        eigv = [eigv EigInfo{i+1}(1)];
    end
end
% semilogy(1:length(plt),plt,'-','linewidth',2)
% % xticks(0:25:length(plt))
% xlabel('Num of QR Factorizations','Interpreter','Latex','FontSize',16)
% ylabel('Residual, $| t_{m,m-1} |$','Interpreter','Latex','FontSize',16)
% hold off

[v,d] = eig(A0);

fprintf('The inf-norm between Matlab ''eig'' command and...\n')

fprintf('\n  ...calculated eigenvalues are:\n')
eigv = sort(eigv)';
for i=1:m
    fprintf('   %.12e   |   %.6e\n',eigv(i),norm(eigv(i)-d(i,i),inf))
end

[ eigvec ] = eigvector( eigv',A0,m );
fprintf('\n\n  ...calculated eigenvectors are:   \n')
for i=1:m
    fprintf('    %.6e    ',norm(eigvec(:,i)-v(:,i),inf))
end
fprintf('\n')
disp(eigvec)

%% Various useful debugging print statements
% [ H, Q ] = tridiag( A0 );
% fprintf('Original Matrix, A0 = \n')
% disp(A0)
% fprintf('Hessenberg Form of A0\n')
% disp(H)
% fprintf('\nunitary Q matrix\n')
% disp(Q)
% fprintf('H = Q''*A0*Q \n')
% disp(Q'*A0*Q)
% fprintf('A0 = Q*dumA*Q'' \n')
% disp(Q*H*Q')

% [ T_new ] = qralg( H );
% fprintf('\nDiag Matrix of Eigenvalues of A0\n')
% disp(T_new) % outputted eigenvalues in diagonal matrix form

