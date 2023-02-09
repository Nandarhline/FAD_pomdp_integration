clear
tic
%%%%%%%%%%%%%%%% Transition matrix definition %%%%%%%%%%%%%
%%%%%%%%%%%% Pablo G. Morato / N. Hlaing %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% 02/03/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Input file: Damage size from MCS
load('dr_FAD.mat')

%%% Input parameters
sta = 40; % Number of damage states
stf = 30;
acrit = max(max(a)); % Critical damage size
samp = size(a,1);
years = size(a,2);

%%% Discretization scheme: definition of the intervals
a0_interv = [1e-20, linspace(0.1603,acrit,sta-1), 1e20];
gf_interv = [-100, linspace(0,2,stf-1), 100];

%%% Counting deterioration rates
%%%Initial state
A0 = a(:,1);
F0 = gFM(:,1);
hist = histcounts2(F0, A0, gf_interv, a0_interv,'Normalization','probability'); % (staq, sta0)
AF0 = reshape(hist, 1, sta*stf); % AQ joint distribution at year 0(f is inner loop)
b0 = zeros(1,sta*stf*years); % Initial belief
b0(1:sta*stf) = AF0; 

%%%Transition probability [T]
ii = 1; %%% Index
T0 = zeros(years*sta*stf,years*sta*stf); %%% Initializing T
for i=1:years-1 % Det. rate index
    A = a(:,i); % Samples a at det. rate i
    A_ = a(:,i+1); % Samples a at det. rate i+1
    F = gFM(:,i); % Samples q at det. rate i
    F_ = gFM(:,i+1); % Samples q at det. rate i+1
    for j=1:sta % a index
        counta = A>=a0_interv(j) & A<a0_interv(j+1);      
        for k=1:stf
            countf = F>=gf_interv(k) & F<gf_interv(k+1);            
            Anext = A_(counta & countf);
            Fnext = F_(counta & countf);
            if sum(counta & countf)<1 % In case there are no samples in that state
                    if sum(countf)==0 && sum(counta)~=0
                        hist = histcounts(A_(counta), a0_interv, 'Normalization','probability');
                        T0(ii,i*sta*stf+1:(i+1)*sta*stf) = repelem(hist,stf)/stf;
                    elseif sum(counta)==0 && sum(countf)~=0
                        hist = histcounts(F_(countf), gf_interv, 'Normalization','probability');
                        T0(ii,i*sta*stf+1:(i+1)*sta*stf) = repmat(hist,1,sta)/sta;
                    else   
                         T0(ii,sta*stf*i+1) = 1;
                    
                    end
            else
            hist = histcounts2(Fnext, Anext, gf_interv, a0_interv, 'Normalization','probability');
            T0(ii,i*sta*stf+1:(i+1)*sta*stf) = reshape(hist, 1, sta*stf);
            end
            
            ii = ii+1;
         end
    end
    
end

%%% Defining the last deterioration rate 
%%% Transition at last det. rate behaves like previous transition
T0((years-1)*sta*stf+1:years*sta*stf,(years-1)*sta*stf+1:years*sta*stf) = ...
    T0((years-2)*sta*stf+1:(years-1)*sta*stf,(years-1)*sta*stf+1:years*sta*stf);
toc

%%% Perfect repair matrix Tr = [b0;b0;...;b0]
Tr = repmat(b0,sta*stf*years,1);
toc

dr_env.T0 = sparse(T0);
dr_env.Tr = Tr;
dr_env.b0 = b0;
dr_env.aint = a0_interv;
dr_env.fint = gf_interv;
dr_env.sta = sta;
dr_env.stf = stf;
save('dr_ENV','dr_env','-v7.3')