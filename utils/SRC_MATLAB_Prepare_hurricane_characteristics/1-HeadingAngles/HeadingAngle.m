clear
clc
close

sample = xlsread('sample.xlsx','Greater Storm');

nsample= length(sample);
ntheta = 73;
omnigreater = 2.3781e-04;
omnilesser  = 2.0647e-04;
hd     = 500;
ho     = 30;

theta  = zeros(ntheta,1);
lamuda = zeros(ntheta,1);

T = max(sample(:,1)) - min(sample(:,1)) + 1;

for i = 1:ntheta
    temp = 0;
    theta(i) = -180+(i-1)*360/(ntheta-1);
    for j = 1:nsample
        temp = temp + exp(-1/2*(sample(j,2)/hd)^2)/sqrt(2*pi)/hd * ...
               exp(-1/2*((sample(j,3)-theta(i))/ho)^2)/sqrt(2*pi)/ho;
    end
    lamuda(i) = temp / T;
end
r=14.269;
t=12.747;
headPrb  = lamuda ./ omnigreater ;
headBeta = ((theta+180)./360).^(r-1).*(1-(theta+180)./360).^(t-1)*gamma(r+t)/gamma(r)/gamma(t)/360/2.15;

% mean = 10.5;
% sd   = 32.6;
% headPrb  = lamuda ./ omnilesser ;
% headnorm = normpdf(theta,mean,sd)/2;



figure(1)
plot(theta,headPrb,'r-o',theta,headBeta,'b--*');grid on;
xlabel('Heading(degrees)')
ylabel('Probabilistic Distribution of Heading Angle')
title('Probabilistic Distribution of Heading Angle for the Greater storms(dP>48mb)')
legend('Calculated Probability','Beta Approximation')

% figure(1)
% plot(theta,headPrb,'r-o',theta,headnorm,'b--*');grid on;
% xlabel('Heading(degrees)')
% ylabel('Probabilistic Distribution of Heading Angle')
% title('Probabilistic Distribution of Heading Angle for the Lesser storms(dP<48mb)')
% legend('Calculated Probability','Normal Approximation')



