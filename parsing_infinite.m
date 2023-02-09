%clear
load('dr_ENV.mat')
T0 = dr_env.T0 ;
Tr = dr_env.Tr;
state0 = dr_env.b0;
sta = dr_env.sta;
stf = dr_env.stf;
PO_D = dr_env.PO_D;
PO_ND = dr_env.PO_ND;

%Costs (failure, repair, inspection)
discRate = 0.94;
cF = -1e6*discRate;
cR = -1.2e4*discRate;
cI = -1e3*discRate*discRate; %Observation is taken later

%%%%%%%%%% Infinite Horizon %%%%%%%%%%
staT = size(T0,1);
rF = zeros(1,staT); % Reward of the state
rF(1:stf:staT) = cF;
Rtrans = rF*T0'; 
R0 = Rtrans - rF; % Immediate reward for action 0 (DN-NI)
R1 = R0 + cI; % Immediate reward for action 1 (DN-I)
R2 = ones(1,staT).*cR; % Immediate reward for action 2 (RP)

dr_env.discRate = discRate;
dr_env.R0 = R0; % Immediate reward for action 0 (DN-NI)
dr_env.R1 = R1; % Immediate reward for action 1 (DN-I)
dr_env.R2 = R2; % Immediate reward for action 2 (RP)
save('IFH_input.mat','dr_env','-v7.3')

