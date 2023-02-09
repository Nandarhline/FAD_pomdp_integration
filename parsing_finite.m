
%%%%%%%%%%%%%%%%%%% POMDP Finite horizon %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Transition matrices %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Pablo G. Morato %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% 09/03/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Loading input file from infinite horizon POMDP
clear
load('IFH_input.mat')
T0 = dr_env.T0;
Tr = dr_env.Tr;
state0 = dr_env.b0;
PO_D = dr_env.PO_D;
PO_ND = dr_env.PO_ND;

%Costs (failure, repair, inspection)
discRate = 0.94;
cF = -1e6*discRate;
cR = -1.2e4*discRate;
cI = -1e3*discRate*discRate; %Observation is taken later

% Number of states for each model
staT = size(T0,1);
nstate = dr_env.sta*dr_env.stf;
yearsT = staT/nstate;

%%% Definition of the time horizon
years = yearsT; %It should be equal to yearsT
%years = 5; % In case a shorter lifetime is desired

%%% Total number of states |St| = |N/2*(S+S*N)| 
%%%(Sum of an arithmetic progression)
expsize = years*(nstate+nstate*years)/2;

% Perfect repair matrix
Sta0Exp = repmat(state0(1:nstate),nstate,1); 

%%% Initial belief
b0fini = zeros(1,expsize);
b0fini(1:nstate) = state0(1:nstate); 
dr_env.b0 = b0fini;

%%% Obervation model %%
loopSize = expsize/nstate; % Same as sum(1:years)
dr_env.PO_D = repmat(PO_D(1:nstate),loopSize,1);
dr_env.PO_ND = repmat(PO_ND(1:nstate),loopSize,1);

%%% Initializing the trainsition matrices 
T0fini = sparse(expsize,expsize); %Initializing T0
Trfini = sparse(expsize,expsize); %Initializing TR

count1 = 1; % Initial index
countn = 1; % Element index initialization
COUNTn = 0; % Row initial index r(t=0) 
COUNTnPlus = 1; % Column initial index c(t=0)

%%% Iterating over the years
for i=1:years-1 %(One year less - Last year wil be defined with an absorbing state)
   innLoops = i; %Inner loops each year
   
   %%% Constructing the transition and reward matrices 
   %%% Iteration depending on the year (and thus deterioration rates)
   for j=1:innLoops
       
       T0fini( (COUNTn+(j-1))*nstate+1:(COUNTn+j)*nstate, (COUNTnPlus+(j))*nstate+1:(COUNTnPlus+j+1)*nstate )...
           = T0( (j-1)*nstate+1:j*nstate, (j)*nstate+1:(j+1)*nstate );

       Trfini( (COUNTn+(j-1))*nstate+1:(COUNTn+j)*nstate, (COUNTnPlus)*nstate+1:(COUNTnPlus+1)*nstate )...
           = Sta0Exp(1:nstate,1:nstate);
   
   end
   
   % Row indices
   countn = count1+(i-1); % r(t) = r(t-1) + 1
   COUNTn = i*((count1+countn)/2); % sum r(t) = t( (r(0) + r(t))/2 )
   
   % Column indices
   countnPlus = count1+(i); % c(t) = c(t-1) + 1
   COUNTnPlus = (i+1)*((count1+countnPlus)/2); % sum c(t) = t( (c(0) + c(t))/2 )
    
end

T0fini( (COUNTn)*nstate+1:(COUNTn+j+1)*nstate,  expsize) = 1; %Absorbing state
Trfini( (COUNTn)*nstate+1:(COUNTn+j+1)*nstate,  expsize) = 1; %Absorbing state

%spy(T0fini); hold on; spy(Trfini)
dr_env.T0 = T0fini;
dr_env.Tr = Trfini;

% Reward for each state
rF = zeros(1,expsize);
rF( 1, 1:dr_env.stf:expsize) =  cF;  
Rtrans = rF*T0fini'; 
R0 = Rtrans - rF;

count1absor = count1+(years-2);
COUNTn = (years-1)*((count1+count1absor)/2);
R0(COUNTn*nstate+1:expsize) = 0;

R1 = R0 + cI;
R1(COUNTn*nstate+1:expsize) = 0;

R2 = ones(1,expsize).*cR;
R2(COUNTn*nstate+1:expsize) = 0;

dr_env.R0 = R0;
dr_env.R1 = R1;
dr_env.R2 = R2;
save('FH_input','dr_env','-v7.3')

