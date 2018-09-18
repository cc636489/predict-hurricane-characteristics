clear
clc
close

sample = xlsread('sample.xlsx','Greater Storm');

%% Initializting
nsample  = length(sample);
ntheta   = 1;
nho      = 29;  
ntry     = 500; 
p        = 0.1; 
hd       = 500;
initial1 = 6; % for ho  
step1    = 6;
initial2 = 0;% for theta
step2    = 6;


CVSErr   = zeros(nho,ntheta); % Directional rate at theta=0;
ho       = zeros(nho,     1);   % Kernel Width ho
theta    = zeros(ntheta,  1);

%% Calculating CVSError
for i = 1:ntheta
    theta(i) = initial2+(i-1)*step2;
    for j = 1:nho
        ho(j) = initial1+(j-1)*step1;
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
                    temp1 = temp1 + exp(-1/2*(sample(k,2)/hd)^2)/sqrt(2*pi)/hd * ...
                        exp(-1/2*((sample(k,3)-theta(i))/ho(j))^2)/sqrt(2*pi)/ho(j);
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
% p = polyfit(ho,CVSErr,4);
% poly2 = polyval(p,ho);

%% FIGURES
figure(1)
for i = 1:ntheta
    plot(ho,CVSErr(:,i),'b--*');hold on;grid on;
end
xlabel('Kernel Size hd(km)')
ylabel('CVSError')
title('CVSE method for the directional storm rate of Greater Storm: theta = 0')

% figure(2)
% plot(hd,poly2(:,1),'b-o',hd,poly2(:,2),'c-*',hd,poly2(:,3),'y-d',hd,poly2(:,4),'r-^',hd,poly2(:,5),'m-s',hd,poly2(:,6),'k-x',hd,poly2(:,7),'g-p',hd,poly2(:,8),'g-<')
% xlabel('Kernel Size hd(km)')
% ylabel('CVSError')
% title('CVSE method for the directional storm rate of Lesser Storm: theta = 0')
% legend('h(alpha)= 12 deg','h(alpha)= 18 deg','h(alpha)= 24 deg','h(alpha)= 30 deg','h(alpha)= 36 deg','h(alpha)= 42 deg','h(alpha)= 48 deg','h(alpha)= 54 deg')




