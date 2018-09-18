clear
clc
close

sample = xlsread('data.xlsx');

%% Initializting
nsample = length(sample); % number of observation
nhd     = 51;  % number of CVSError point on one line
ntry    = 500; % number of partitions
p       = 0.1; % Probability for validation samples
initial = 100; % Initial of hd
step    = 20;  % step of hd

CVSErr  = zeros(nhd,1);   % Omni-directional rate
hd      = zeros(nhd,1);   % Kernel Width hd

%% Calculating CVSError
for i = 1:nhd
    hd(i) = initial+(i-1)*step;
    temperror = 0;
    for j = 1:ntry
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
                temp1 = temp1 + exp(-1/2*(sample(k,2)/hd(i))^2)/sqrt(2*pi)/hd(i);
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
    CVSErr(i) = temperror;
end

%% FIGURES
figure(1)
plot(hd,CVSErr,'b--o');grid on;
xlabel('Kernel Size hd(km)')
ylabel('CVSError')
title('CVSE method for the omni-directional storm rate of Lesser Storm')






