clear
clc
close
tic

%% Initializing
global resample nresample hd

sample = xlsread('sample.xlsx','Greater Storm');
nsample = length(sample);

hd       = 500;
nboot    = 1000;
poismean = 1;

%% Calculating u & k 
utemp = 0;
ktemp = 0;
for i = 1:nboot
    poisnum = poissrnd(poismean,1,nsample);%1 means 1row
    nresample = sum(poisnum);
    resample  = zeros(nresample,2);
    %%%%%%%% Constructing a bootstrapping resample set
    temp = 0;
    for j = 1:nsample
        if ( poisnum(j) ~= 0 )
            k = 1;
            while ( k <= poisnum(j) )
                resample(k+temp,1) = sample(j,2);%distance
                resample(k+temp,2) = sample(j,4);%distance
                k = k+1;
            end
            temp = temp + poisnum(j);
        end
    end
    %%%%%%%% Optimization to obtain u and k
    [ubar,kbar] = optimukfun();
    utemp = utemp + ubar;
    ktemp = ktemp + kbar;
    disp(['Finish bootstap ',num2str(i),' already !']);
end

uaverage = utemp / nboot;
kaverage = ktemp / nboot;

uaverage
kaverage

toc
