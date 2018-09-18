clear
clc
close

sample1 = xlsread('sample.xlsx','Greater Storm');
nsample1 = length(sample1);

%pressure=sample(:,4)
%radius=sample(:,5)

increase1 = sort(sample1(:,5));
decrease1 = sort(sample1(:,5),'descend');
per16cent1 = increase1( ceil(nsample1*0.16) );
per84cent1 = decrease1( ceil(nsample1*0.16) );

data16set1 = zeros(ceil(nsample1*0.16),2);
data84set1 = zeros(ceil(nsample1*0.16),2);
medianset1 = zeros(nsample1-ceil(nsample1*0.16)*2,2);

k = 1;
m = 1;
j = 1;
for i = 1:nsample1
    if ( sample1(i,5) <= per16cent1 && k <= ceil(nsample1*0.16))
        data16set1(k,1)=sample1(i,4);
        data16set1(k,2)=sample1(i,5);
        k = k+1;
    elseif ( sample1(i,5) >= per84cent1 && m <= ceil(nsample1*0.16))
        data84set1(m,1)=sample1(i,4);
        data84set1(m,2)=sample1(i,5);
        m = m+1;
    else
        medianset1(j,1)=sample1(i,4);
        medianset1(j,2)=sample1(i,5);
        j = j+1;
    end 
end

logdata16set1 = log(data16set1);
logdata84set1 = log(data84set1);
logmedianset1 = log(medianset1);
logdpdataset1 = log(sample1(:,4));
logrpdataset1 = log(sample1(:,5));

p16coeff1 = polyfit(logdata16set1(:,1),logdata16set1(:,2),1);
p84coeff1 = polyfit(logdata84set1(:,1),logdata84set1(:,2),1);
palcoeff1 = polyfit(logdpdataset1,logrpdataset1,1);

figure(1)
dp16val1 = 40:1:120;
rp16val1 = exp( polyval(p16coeff1,log(dp16val1)) );
rp84val1 = exp( polyval(p84coeff1,log(dp16val1)) );
rpmeanval1   = (rp16val1 + rp84val1)/2;
mediancoeff  = polyfit(log(dp16val1),log(rpmeanval1),1);
rpmedianval1 = exp( polyval(palcoeff1,log(dp16val1)) );

plot(dp16val1,rp16val1);hold on;
plot(dp16val1,rpmeanval1);hold on;
plot(dp16val1,rpmedianval1);hold on;
plot(dp16val1,rp84val1);hold on;
scatter(sample1(:,4),sample1(:,5),'filled','b');grid on;

title('Model:ln(Rp)= 3.4887-0.0847ln(Dp);Greater Storm');
xlabel('dp(mb)');
ylabel('Rp(nmile)');
legend('16%','mean','median','84%','Original Data');



