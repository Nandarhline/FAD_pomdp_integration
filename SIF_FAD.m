% Probabilistic Analysis using SIF Model with Failure Assessment Diagram % %
% Nandar Hlaing
clear 
tic

% MCS 1e4 samples 60 seconds
% MCS 1e5 samples 700 seconds
% MCS 1e6 samples 17950 seconds

%%% Input Data %%%
aa = [0.1603 -27.6302 0.4599]; % or load the calibration file load('Values_SIF.mat')
n = 3.5e7; % number of cycles per year
T = 20; % lifetime
samp = 1e5;
q_mean = 6.4839; % weibull scale parameter
q_cov = 0.2;
h = 0.8; % weibull shape parameter
DOB = 0.81;% Degree of bending (Stress Ratio for bending)
acrit = 16; % critical crack size (Thickness, mm) 

%% Respective Distributions of Uncertainties %%
% Initial crack size as exponential distribution %
a0_mean = aa(1); % initial crack size (mm)
a0 = exprnd(a0_mean,samp,1);
r = 0.2;
c0 = a0./r;

% Paris law parameters as normal distributions %
m = 3;
lnCa_mean = aa(2);
lnCa_std = aa(3);
Ca = lognrnd(lnCa_mean,lnCa_std,samp,1);

% Stress range as normal distribution %
q_std = q_cov*q_mean;
q = normrnd(q_mean,q_std,samp,1);
S = q.*gamma(1+1/h);

% Taking out negative samples 
S_check=mean(S(S<0));
while S_check<0
    S(S<0)=normrnd(q_mean,q_std,length(S_check),1);
    S_check=mean(S(S<0));
end

% Maximum primary stress as normal distribution %
S_mean = q_mean.*gamma(1+1/h);
Ps = S/2;
Pb = DOB.*Ps; % Bending stress
Pm = (1-DOB).*Ps; % Membrane stress

% Yield strength as lognormal distribuiton %
yield_median = 355; % Yield strength
yield_cov = 0.07;
yield_mean = log(yield_median); % Yield strength
yield_std = (log(yield_cov^2+1))^0.5;
yield = lognrnd(yield_mean,yield_std,samp,1);

% Fracture toughness as 3 parameters weibull distrition %  
Kmat_k = 4; % shape parameter
Kmat_theta = 20; % threshold parameter
T1 = 10; % to calculate scale parameter
T_27J = 20;
Kmat_sigma = (11+77*exp((T1-T_27J-7)/52))*(25/acrit)^0.25*(log(1/0.95))^0.25; 
F = wblrnd(Kmat_sigma,Kmat_k,samp,1);
Kmat = (F+Kmat_theta);
kmean = mean(Kmat);

% Resistance fracture toughness as lognormal distribution %
Rf_median = 1.7;
Rf_cov = 0.18;
Rf_mean = log(Rf_median); % Yield strength
Rf_std = (log(Rf_cov^2+1))^0.5;
Rf = lognrnd(Rf_mean,Rf_std,samp,1);

% Residual stress as lognormal distribution %
Rs_median = 300;
Rs_cov = 0.2;
Rs_mean = log(Rs_median); % Yield strength
Rs_std = (log(Rs_cov^2+1))^0.5;
Rs = lognrnd(Rs_mean,Rs_std,samp,1);

%% Calculating Crack Growth %%
a = zeros(samp,T);
c = zeros(samp,T);
% Lr = zeros(samp,T);
% Kr = zeros(samp,T);
% rho = zeros(samp,T);
Mkma = zeros(samp,T+1);
Mkba = zeros(samp,T+1);
Mkmc = zeros(samp,T+1);
Mkbc = zeros(samp,T+1);
Yma = zeros(samp,T+1);
Yba = zeros(samp,T+1);
Ymc = zeros(samp,T+1);
Ybc = zeros(samp,T+1);


parfor j=1:samp
    Yinit = [a0(j) c0(j)];
    Ca1 = Ca(j);
    Cc1 = Ca(j);
    S1 = S(j);
    acheck = a0(j);
    %[Mkmat,Mkbat,Mkmct,Mkbct,Ymat,Ybat,Ymct,Ybct] = geom_magtub(Yinit(1),Yinit(2),acrit,acheck);
    for i=1:T
        N = [1+(i-1)*n i*n];
        [N,Y] = ode45(@(N,Y) CG_Func(N,Y,Ca1,Cc1,m,DOB,S1,acrit,acheck),N,Yinit);
        %[N,Y] = ode45(@(N,Y) CG_Func1(N,Y,Ca1,Cc1,m,DOB,S1,Mkmat,Mkbat,Mkmct,Mkbct,Ymat,Ybat,Ymct,Ybct),N,Yinit);
        
        a(j,i) = abs(Y(end,1));
        c(j,i) = abs(Y(end,2));
        
         % Remove instabilites for failed crack
        if  a(j,i)<Yinit(1) || a(j,i)>acrit || isnan(a(j,i)) 
            a(j,i) = acrit+0.1;
        end
        
        if  c(j,i)<Yinit(2) || isnan(c(j,i))
            c(j,i) = abs(max(Y(:,2))); 
        end
        
        if c(j,i)>160
            c(j,i)=160;
        end
        
    Yinit = [a(j,i) c(j,i)];
    end    
end

a = [a0 a];
c = [c0 c];

for j=1:samp
[Mkma(j,:),Mkba(j,:),Mkmc(j,:),Mkbc(j,:),Yma(j,:),Yba(j,:),Ymc(j,:),Ybc(j,:)]=geom_magtub(a(j,:),c(j,:),acrit,a0(j));
end

%% Calculating failure probability
for i=1:T+1
    Ki_s = abs(Rs.*(Yma(:,i).*Mkma(:,i).*(1-DOB) + Yba(:,i).*Mkba(:,i).*DOB).*sqrt(pi.*a(:,i)*0.001));
    Ki_p = abs(Ps.*(Yma(:,i).*Mkma(:,i).*(1-DOB) + Yba(:,i).*Mkba(:,i).*DOB).*sqrt(pi.*a(:,i)*0.001));
%       Ki_s = Rs.*sqrt(pi.*a(:,i).*0.001);
%       Ki_p = Ps.*sqrt(pi.*a(:,i).*0.001);
     Ki = Ki_p+Ki_s;
     
     alpha = (a(:,i)./acrit)./(1+(acrit./c(:,i)));
     sigma_ref = (Pb+((Pb.^2 + 9.*Pm.^2.*(1-alpha).^2)).^0.5)./(3.*(1-alpha).^2);
     Lrt = sigma_ref./yield;
    
     zeta = Ki_s.*Lrt./Ki_p;
     rho1 = 0.1*zeta.^0.714-0.007*zeta.^2+3e-5*zeta.^2;
     rho1(zeta>4) = 0.25;
    
    rho = 4.*rho1.*(1.05-Lrt);
    rho(Lrt<=0.8) = rho1(Lrt<=0.8);
    rho(Lrt>=1.05) = 0;
    
    Krt = Ki./Kmat + rho;
    gFM(:,i) = Rf-sqrt(Krt.^2+Lrt.^2);
    pf(i) = sum(gFM(:,i)<0)/samp;
    Kr(:,i) = Krt;
    Lr(:,i) = Lrt;
    rhot(:,i) = rho;

end
beta_FAD = -norminv(pf);

save('dr_FAD.mat','a','gFM')
save('Crack_FAD.mat','a','c','q','Ca','yield','Kmat','Rs','Rf','Kr','Lr','gFM','rhot')
save('pf_FAD.mat','pf')

toc



