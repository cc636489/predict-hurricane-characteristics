function [c,ceq] = mycon(x0)
% u: x0(1,1)
% k: x0(2,1)
c = real( x0(1,1)*(x0(2,1)-1)/x0(2,1) - (48/x0(1,1))^x0(2,1));
ceq=[];

end