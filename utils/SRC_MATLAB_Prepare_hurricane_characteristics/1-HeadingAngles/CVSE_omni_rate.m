clear
clc
close

sample    = xlsread('sample.xlsx','Lesser Storm');
nsample   = length(sample); % number of observation
omnitemp  = 0;
hd        = 500;
max       = sample(1,1);
min       = sample(1,1);

%% Calculating omni-directional rate

for k = 1:nsample
    omnitemp = omnitemp + exp(-1/2*(sample(k,2)/hd)^2)/sqrt(2*pi)/hd;
    if (sample(k,1) > max)
        max = sample(k,1);
    elseif (sample(k,1) < min)
        min = sample(k,1);
    end
end
T    = max - min + 1;
omni = omnitemp / T;

disp(['omni-directional rate is',num2str(omni),'storms/km/year']);
