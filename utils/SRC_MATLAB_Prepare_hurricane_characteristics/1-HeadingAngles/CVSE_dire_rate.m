clear
clc
close

sample = xlsread('sample.xlsx','Lesser Storm');

nsample= length(sample);
ntheta = 73;
hd     = 500;
ho     = 30;

lamuda = zeros(ntheta,1);
theta  = zeros(ntheta,1);

T = max(sample(:,1)) - min(sample(:,1)) + 1;

for i = 1:ntheta
    temp = 0;
    theta(i) = -180+(i-1)*360/(ntheta-1);
    for j = 1:nsample
        temp = temp + exp(-1/2*(sample(j,2)/hd)^2)/sqrt(2*pi)/hd * ...
               exp(-1/2*( ( sample(j,3)-theta(i) ) / ho)^2)/sqrt(2*pi)/ho;
    end
    lamuda(i) = temp / T;
end

figure(1)
plot(theta,lamuda,'r-o');grid on;
xlabel('Heading(degrees)')
ylabel('Directional Rate (storms/km/deg/year)')
title('Directional rates of storm heading for the Lesser storms(dP<48mb)')


