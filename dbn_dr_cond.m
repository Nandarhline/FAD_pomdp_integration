clear

%%% Input file: DBN tuple
load('dr_ENV.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nstates = length(dr_env.b0);
years = nstates/dr_env.sta/dr_env.stf;
sta = dr_env.sta;
stf = dr_env.stf;

% %% Observation Matrix for Inspection
aPOD = (dr_env.aint(1:end-1)+dr_env.aint(2:end))./2;
aPOD(end) = dr_env.aint(end-1)+1;
% %%% DNV-GL inspection model
B0 = 7.3074;
B1 = 2.092; 
eta_sig = 4.189;
thres = 5.4898;

pod = 1-normcdf((thres-(B0+B1*log(aPOD)))/eta_sig);  
POD_pdf1 = repelem(pod,1,stf);
dr_env.PO_D = repmat(POD_pdf1',years,1);
dr_env.PO_ND = 1 - dr_env.PO_D;

% X0 = 0.45; % Normal conditions 
% b = 0.9;
% pod = 1 - 1./ (1 + (aPOD./X0).^b );
% POD_pdf1 = repelem(pod,1,stq);
% dr_env.PO_D = repmat(POD_pdf1',lifet,1);
% dr_env.PO_ND = 1 - dr_env.PO_D;
save('dr_ENV','dr_env','-v7.3')

%%% Prescribing observations
Ins_int = [8]; % Specific observation
Ins_int = Ins_int+1;
% Ins_inter = 5; % observation interval
% Ins_int = Ins_inter+1:Ins_inter:lifet+1;

ins = zeros(1,years+1); % Initializing observation vector
ins(Ins_int) = 1; % Setting the observations

pF = zeros(1,years);
pF(1) = sum(dr_env.b0(1:stf:sta*stf));
pF_plot = zeros(1,years+sum(ins)); % PF for plotting
t_plot = zeros(1,years+sum(ins)); % t for plotting
t_ind = 2; % Index for plotting
b0 = dr_env.b0;
for t=2:years
   state = b0*dr_env.T0;
   AF = reshape(state((t-1)*(sta*stf)+1:t*(sta*stf)),stf,sta); %(stq,sta)
   pF(t) = sum(AF(1, :));
   t_plot(t_ind) = t_plot(t_ind-1) + 1; % t for plotting
   pF_plot(t_ind) = state(t*dr_env.sta); % PF for plotting
   
      if ins(t)==1
       
       state = state.*dr_env.PO_ND';
       normaliz = sum(state);
       state = state./normaliz;

       AF = reshape(state((t-1)*(sta*stf)+1:t*(sta*stf)),stf,sta); %(stq,sta)
       pF(t) = sum(AF(1, :));
       t_ind = t_ind + 1; % Index for plotting purposes
       t_plot(t_ind) = t_plot(t_ind-1); % t for plotting
       pF_plot(t_ind) = state(t*dr_env.sta); % PF for plotting 
      end
   
   b0 = state;
   t_ind = t_ind + 1; % Index for plotting purposes
end

pf_DBN2 = diff(pF); %Annual failure probability
beta_DBN2 = -norminv(pf_DBN2); %Reliability index
pF_DBN2 = pF;

tt = 1:years-1;
%Plotting
plot(tt,beta_DBN2,'r','LineWidth',1.5)
hold all
xlabel('Time (years)')
ylabel('\beta')
xlim([0 20])
grid minor

%%%Comparison with FMs
load('mcs_upd_out.mat')
%  Scaling
pF_est = (pF_DBN2 - mean(PF_MCS)) ./ std(PF_MCS);
pF_FM = (PF_MCS - mean(PF_MCS)) ./ std(PF_MCS);
MSE = sum((pF_FM - pF_est).^2)/length(pF_est);

%%% Uncertainty of the estimator
std_estim = sqrt(PF_MCS.*(1-PF_MCS)./1e7); 
plot(tt,beta_MCS,'b--','LineWidth',1.5)
legend('DBN-NS','MCS','Location','northwest')



