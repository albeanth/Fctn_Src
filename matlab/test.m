clc,clear
m = 32;
fprintf('')
Full = 0;
Hess = 0;
for idx = 1:m-1
    Full = Full+ (m-idx+1)*4;
    Hess = Hess + 2*4;
    
end
% Full
% Hess