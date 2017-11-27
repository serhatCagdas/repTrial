clear all;
close all;
clc;

%% var covariance test
mu = [0 0];
RL = [0.1 0.3 0.5 0.7 0.9];
NL = [30 50 100 250 500 1000];
% NL = [250 500 1000];
% RL = [0.95  0.9 0.95];

TN = 100;   %%Trial num;

Ig = -0.5*log(1-RL.^2);

for i = 1:length(RL)
    for j = 1:length(NL)
        for k = 1:TN 
            r = RL(i);
            N = NL(j);
            Sigma = [1 r; r 1];
            data  = mvnrnd(mu,Sigma,N); 
            T = table(data);
            %writetable(T,'dataset.xlsx','Sheet',i,'Range','A1');
            X{i,j,k} = data;
        end
    end
end

Icd     = zeros(length(RL),length(NL),TN);
Iuv     = zeros(length(RL),length(NL),TN);
Ikr     = zeros(length(RL),length(NL),TN);
time_CF = zeros(length(RL),length(NL));
time_UV = zeros(length(RL),length(NL));
time_KR = zeros(length(RL),length(NL));


for j = 1:length(NL)
    for i = 1:length(RL) 
        tic
        for k = 1:TN  
            Iuv(i,j,k)   = MutualInformation_uvParam(X{i,j,k});
        end
        time_UV(i,j) = toc; 
        
        tic
        for k = 1:TN  
            Icd(i,j,k)   = conditional_dependency4( X{i,j,k} );
        end
        time_CF(i,j) = toc;
        
        tic
        for k = 1:TN  
            Ikr(i,j,k)   = kraskov_MI_light(X{i,j,k}); 
        end
        time_KR(i,j) = toc; 
    end 
% NL(j)
end



err_CF = zeros(length(RL),length(NL));
std_CF = zeros(length(RL),length(NL));
err_UV = zeros(length(RL),length(NL));
std_UV = zeros(length(RL),length(NL));
err_KR = zeros(length(RL),length(NL));
std_KR = zeros(length(RL),length(NL));

for i = 1:length(RL)
    Icdr=reshape(Icd(i,:,:),length(NL),TN);
    Iuvr=reshape(Iuv(i,:,:),length(NL),TN);
    Ikrr=reshape(Ikr(i,:,:),length(NL),TN);
    err_CF(i,:) = mean((Icdr' - Ig(i)).^2);
    std_CF(i,:) = std((Icdr' - Ig(i)).^2);
    err_UV(i,:) = mean((Iuvr' - Ig(i)).^2);
    std_UV(i,:) = std((Iuvr' - Ig(i)).^2);
    err_KR(i,:) = mean((Ikrr' - Ig(i)).^2);
    std_KR(i,:) = std((Ikrr' - Ig(i)).^2);
end
