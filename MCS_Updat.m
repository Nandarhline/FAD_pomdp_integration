clear
%%%%%%%%%%%%%%%%%%% Reliability updating %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Pablo G. Morato / N. Hlaing %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% 02/03/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Input file: MCS 
load('dr_FAD.mat')
[samp, years] = size(a); % Importing lifet and samples
%acrit = max(max(aa))-0.1; % Critical damage size %% To be removed?

%%% Prescribing observations
Ins_int = [8 ]; % Specific observation
Ins_int = Ins_int+1;

% Ins_inter = 5; % observation interval
% Ins_int = Ins_inter+1:Ins_inter:lifet+1;

ins = zeros(1,years+1); % Initializing observation vector
ins(Ins_int) = 1; % Setting the observations

%%% Discretization of inspection data
%%% DNV-GL inspection model
B0 = 7.3074;
B1 = 2.092; 
eta_sig = 4.189;
thres = 5.4898;


%%% COMPUTATION
PF_MCS = zeros(1,years); % Initializing PF
PoD_ = ones(samp,1); % Initializing observation conditionals
pF_plot = zeros(1,years+sum(ins)); % PF for plotting
t_plot = zeros(1,years+sum(ins)); % t for plotting
t_ind = 1; % Index for plotting
for t=1:years
    
    at = a(:,t); % a(t) 
    gf = gFM(:,t) <= 0; % Fail conditional: a>ac   
    Num = sum(gf.*PoD_)/samp; % P(a<ac & a>ad)
    Den = sum(PoD_)/samp; % P(a>ad) 
    PF_MCS(1,t) = Num/Den; % P(a<ac & a>ad)/P(a>ad)    
    pF_plot(t_ind) = Num/Den; % PF for plotting
    t_plot(t_ind) = t-1; % t for plotting

    if ins(t) == 1 
        %%% DNV-GL observation model
        UnifRand = rand(samp,1);
        pod = 1-normcdf((thres-(B0+B1*log(a(:,t))))/eta_sig);   
        gH = UnifRand >= pod; % Observation conditional: u>=pod for not detected
        PoD_ = PoD_ .* gH; % Combining observation conditions        
        Num = sum(gf.*PoD_)/samp; % P(a<ac & a>ad)
        Den = sum(PoD_)/samp; % P(a>ad)
        PF_MCS(1,t) = Num/Den; % P(a<ac & a>ad)/P(a>ad) 
        
        t_ind = t_ind + 1; % Index for plotting purposes
        t_plot(t_ind) = t-1; % t for plotting
        pF_plot(t_ind) = Num/Den; % PF for plotting 
    end
    
    t_ind = t_ind + 1; % Index for plotting purposes    
end

%%% Post-processing
tt = 1:years-1; % Time steps
pf_MCS = diff(PF_MCS); % Computation of Pf (Annual failure probability)
beta_MCS = -norminv(pf_MCS); %Reliability index
std_estimPF = sqrt(PF_MCS.*(1-PF_MCS)./samp); % Uncertainty of the estimator PF
std_estimpf = sqrt(pf_MCS.*(1-pf_MCS)./samp); % Uncertainty of the estimator Pf
plot(tt,beta_MCS,'r','LineWidth',1.5)
hold all

%%% Saving file
save('mcs_upd_out.mat','PF_MCS','pf_MCS','beta_MCS') % Saving reliability