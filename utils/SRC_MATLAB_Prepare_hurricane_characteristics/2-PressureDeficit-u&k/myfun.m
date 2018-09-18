function [myval,temp] = myfun(x0)% x is a parameter vector(u,k), r is barrier parameter.

% x0 = [1000;5];
global resample nresample hd
% resample  = xlsread('sample.xlsx','Greater Storm');
% nresample = length(resample);
% hd        = 500;
% x0   = [81.4724;4.5290];
temp = 0;
for i = 1:nresample
    temp = temp + exp( -(resample(i,1)/hd)^2/2 )/sqrt(2*pi)/hd  * log( (x0(2,1)/x0(1,1))* ...
          (resample(i,2)/x0(1,1))^(x0(2,1)-1)*exp( -(resample(i,2)/x0(1,1))^x0(2,1)+(48/x0(1,1))^x0(2,1)));
end
myval = - real(temp);

end
