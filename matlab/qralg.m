function [ D_new,epsi ] = qralg( T,flag )
% Take triangular (Hessenberg) matrix and find eigenvalue(s) of A using
% either unshifted QR (flag=0) or Wilkinson Shifted QR (flag=1).
% Also returns the residual as a function of QR iterations required for  
% convergence.

D_new = T;
m = length(D_new);
epsi = [1];
cnt = 0;
while epsi(end) >= 1E-12
   cnt = cnt + 1;
   D_old = D_new;
   
   if flag == 0 % unshifted QR
       [q,r] = qr(D_old); % QR factorization of hessenberg matrix
       D_new = r*q;
   elseif flag == 1 % Wilkinson shifted QR
       B = D_new(end-1:end,end-1:end);
       del = (B(1,1)-B(2,2))/2;
       if del == 0
           tmpsign = 1;
       else
           tmpsign = sign(del);
       end
       mu = B(2,2)-(tmpsign*B(2,1)^2)/(abs(del)+sqrt(del^2+B(2,1)^2));
       [q,r] = qr(D_old-mu*eye(m));
       D_new = r*q + mu*eye(m);
   end
   
   if epsi(end)==1
       epsi = [abs(D_new(m,m-1))];
   else
       epsi = [epsi abs(D_new(m,m-1))];
   end
end

end