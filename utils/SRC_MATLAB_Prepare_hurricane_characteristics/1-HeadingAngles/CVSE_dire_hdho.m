clear
clc
close

sample = xlsread('sample.xlsx','Lesser Storm');

%% Initializting
nsample  = length(sample); % number of observation
nhd      = 51;  % number of CVSError point on one line
nho      = 8;   % number of ho lines
ntry     = 200; % number of partitions
p        = 0.5; % Probability for validation samples
initial1 = 100; % Initial of hd
initial2 = 12;  % Initial of ho
step1    = 20;  % step of hd
step2    = 6;   % step of ho
theta    = 0;   % Direction of interest

CVSErr   = zeros(nhd,nho); % Directional rate at theta=0;
hd       = zeros(nhd,1);   % Kernel Width hd
ho       = zeros(nho,1);   % Kernel Width ho
poly2    = zeros(nhd,nho);

%% Calculating CVSError
for i = 1:nho
    ho(i) = initial2+(i-1)*step2;
    for j = 1:nhd
        hd(j) = initial1+(j-1)*step1;
        temperror = 0;
        for pp = 1:ntry
            %%%%%%%%%%%% Initializing
            temp1   = 0;
            temp2   = 0;
            max1    = 0;
            max2    = 0;     %%% Intentionally set zero  years in order to go into if
            min1    = 3000;
            min2    = 3000;  %%% Intentionally set large years in order to go into if
            %%%%%%%%%%%% Random partition scheme
            [train, test] = crossvalind('HoldOut', nsample, p);
            ntrain  = length(train);
            ntest   = length(test);
            %%%%%%%%%%%% Building training set
            for k = 1:ntrain
                if(train(k)==0)
                    continue
                else
                    max1 = sample(k,1);
                    min1 = sample(k,1);
                end
                break
            end
            for k = 1:ntrain
                if(train(k)~=0)
                    temp1 = temp1 + exp(-1/2*(sample(k,2)/hd(j))^2)/sqrt(2*pi)/hd(j) * ...
                        exp(-1/2*((sample(k,3)-theta)/ho(i))^2)/sqrt(2*pi)/ho(i);
                    if (sample(k,1) > max1)
                        max1 = sample(k,1);
                    end
                    if (sample(k,1) < min1)
                        min1 = sample(k,1);
                    end
                end
            end
            %%%%%%%%%%%% Building testing set
            for l = 1:ntest
                if(train(l)==0)
                    continue
                else
                    max2 = sample(l,1);
                    min2 = sample(l,1);
                end
                break
            end
            for l = 1:ntest
                if (test(l)~=0 && sample(l,2)<=40)
                    temp2 = temp2 + 1;
                    if (sample(l,1) > max2)
                        max2 = sample(l,1);
                    end
                    if (sample(l,1) < min2)
                        min2 = sample(l,1);
                    end
                end
            end
            
            T1 = max1 - min1 + 1;
            T2 = max2 - min2 + 1;
            
            esti_rate = temp1/T1*(1-p);
            vali_rate = temp2/80/T2*p;
            temperror = temperror + (esti_rate - vali_rate)^2;
        end
        CVSErr(j,i) = temperror;
    end
end

%% Polynomial Regression
for i = 1:nho
    p = polyfit(hd,CVSErr(:,i),4);
    poly2(:,i) = polyval(p,hd);
end

%% FIGURES
figure(1)
plot(hd,poly2(:,1),'b-o',hd,poly2(:,2),'c-*',hd,poly2(:,3),'y-d',hd,poly2(:,4),'r-^',hd,poly2(:,5),'m-s',hd,poly2(:,6),'k-x',hd,poly2(:,7),'g-p',hd,poly2(:,8),'g-<')
xlabel('Kernel Size hd(km)')
ylabel('CVSError')
title('CVSE method for the directional storm rate of Lesser Storm: theta = 0')
legend('h(alpha)= 12 deg','h(alpha)= 18 deg','h(alpha)= 24 deg','h(alpha)= 30 deg','h(alpha)= 36 deg','h(alpha)= 42 deg','h(alpha)= 48 deg','h(alpha)= 54 deg')




