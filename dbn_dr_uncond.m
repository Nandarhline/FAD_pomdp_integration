clear
tic
load('dr_ENV.mat')

sta = dr_env.sta;
stf = dr_env.stf;
b0 = dr_env.b0;
nstates = length(b0);
lifet = nstates/(sta*stf);

pF = zeros(1,lifet);
for t=2:lifet
   state = b0*dr_env.T0;
   AF = reshape(state((t-1)*(sta*stf)+1:t*(sta*stf)),stf,sta); %(stf,sta)
   pF(t) = sum(AF(1, :));
   b0 = state;
end
toc

pf_DBN = diff(pF); %Annual failure probability
%beta_DBN = -norminv(pf_DBN); %Reliability index
pF_DBN = pF;

tt = 0:lifet-1;

%Plotting
semilogy(tt,pF_DBN,'r','LineWidth',1.5)
hold all
% plot(tt,beta_FM,'LineWidth',1.5)
xlabel('Time (years)') 
ylabel('P_F')
grid minor

%%%Comparison with FMs
load('pf_FAD.mat')
semilogy(tt,pf,'b--','LineWidth',1.5)
legend('DBN-DR','MCS','Location','northwest')



