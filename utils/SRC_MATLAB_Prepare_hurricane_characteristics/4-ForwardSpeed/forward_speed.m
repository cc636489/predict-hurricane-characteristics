clear
clc
close

sample1 = xlsread('sample.xlsx','Greater Storm');
sample2 = xlsread('sample.xlsx','Lesser Storm');
%forwardspeed:sample(i,6)

speed  = 1:0.1:16;

mean1   = 1.839069995;% Greater
stdevi1 = 0.437116311;% Greater
lognormal1 = exp(-(log(speed)-mean1).^2./2./stdevi1.^2)./speed./stdevi1./sqrt(2*pi);

mean2   = 1.612184009; % lesser
stdevi2 = 0.435724715; % lesser
lognormal2 = exp(-(log(speed)-mean2).^2./2./stdevi2.^2)./speed./stdevi2./sqrt(2*pi);

figure(1)
h1 = histogram(sample1(:,6),'Normalization','pdf');
h1.BinEdges = [1:16];
hold on;
plot(speed,lognormal1,'r-o');
xlabel('Vf(m/s)');
ylabel('Probability Density');
title('Forward speed probability density function for Greater Storm');


figure(2)
h2 = histogram(sample2(:,6),'Normalization','pdf');
h2.BinEdges = [1:16];
hold on;
plot(speed,lognormal1,'r-o');
xlabel('Vf(m/s)');
ylabel('Probability Density');
title('Forward speed probability density function for Lesser Storm');

