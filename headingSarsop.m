clear


head01 = fopen('../05_Policies/FADfini.pomdp','w'); %opening the file

% Run parsing01 or load the dr_ENV.mat file
%parsing01;

load('FH_input.mat')
T0 = dr_env.T0;
Tr = dr_env.Tr;
state0 = dr_env.b0;
PO_D = dr_env.PO_D;
PO_ND = dr_env.PO_ND;
R0 = dr_env.R0;
R1 = dr_env.R1;
R2 = dr_env.R2;
discRate = dr_env.discRate;

% T0 = full(T0);
% Tr = full(Tr);

[rpa,clpa] = size(T0);

%%% Heading
fprintf(head01,'\n');
fprintf(head01,'%s %f\n','discount:',discRate);
fprintf(head01,'%s %1s\n','values:','reward');
fprintf(head01,'%s %d\n','states:',rpa);
fprintf(head01,'%s %d\n','actions:',3);
fprintf(head01,'%s %d\n','observations:',2);

%%% Initial state definition
fprintf(head01,'\n');
fprintf(head01,'%s','start:');
for i=1:rpa
   %fprintf(head01,' %31.30f',state0(i)); 
   fprintf(head01,' %11.6f',state0(i)); 
end
fprintf(head01,'\n');

tic
%%% Transition matrix Action 0 (Do-nothing - No inspection)
fprintf(head01,'\n');
for i=1:rpa
    T0full = full(T0(i,:));
    for j=1:clpa        
       if T0full(j)~=0 
       fprintf(head01,'%s %-6d %s %-6d %s %-6d %11.6f\n','T:',0,':',i-1,':',j-1,T0full(j)); 
       end
    end
end
toc 

%%% Transition matrix Action 1 (Do-nothing - Inspection)
fprintf(head01,'\n');
for i=1:rpa
    T0full = full(T0(i,:));
    for j=1:clpa        
       if T0full(j)~=0 
       fprintf(head01,'%s %-6d %s %-6d %s %-6d %11.6f\n','T:',1,':',i-1,':',j-1,T0full(j)); 
       end
    end
end
toc

%%% Transition matrix Action 2 (Repair)
fprintf(head01,'\n');
for i=1:rpa
    Trfull = full(Tr(i,:));
    for j=1:clpa        
       if Trfull(j)~=0 
       fprintf(head01,'%s %-6d %s %-6d %s %-6d %11.6f\n','T:',2,':',i-1,':',j-1,Trfull(j)); 
       end
    end
end
toc

%%% Emission matrix Action 0 (No-inspection)
fprintf(head01,'\n');
for i=1:rpa
fprintf(head01,'%s %-6d %s %-6d %s %-6d %11.6f\n','O:',0,':',i-1,':',0,1); 
end

%%% Emission matrix Action 1 (Inspection - PoD)
fprintf(head01,'\n');
for j=1:clpa        
    fprintf(head01,'%s %-6d %s %-6d %s %-6d %11.6f\n','O:',1,':',j-1,':',0,PO_ND(j)); 
end
fprintf(head01,'\n');
for j=1:clpa        
    fprintf(head01,'%s %-6d %s %-6d %s %-6d %11.6f\n','O:',1,':',j-1,':',1,PO_D(j)); 
end

%%% Emission matrix Action 2 (No-inspection)
fprintf(head01,'\n');
for i=1:rpa
fprintf(head01,'%s %-6d %s %-6d %s %-6d %11.6f\n','O:',2,':',i-1,':',0,1); 
end

%%% Reward Action 0 (Do-nothing - No Inspection)
fprintf(head01,'\n');
for i=1:rpa       
       if R0(i)~=0
    fprintf(head01,'%s %-6d %s %-6d %s %s %s %s %-10.1f\n',...
    'R:',0,':',i-1,':','*',':','*',R0(i));
       end
end

%%% Reward Action 1 (Do-nothing - Inspection)
fprintf(head01,'\n');
for i=1:rpa       
       if R1(i)~=0
    fprintf(head01,'%s %-6d %s %-6d %s %s %s %s %10.1f\n',...
    'R:',1,':',i-1,':','*',':','*',R1(i));
       end
end

%%% Reward Action 2 (Repair - No-inspection)
fprintf(head01,'\n');
for i=1:rpa       
       if R2(i)~=0
    fprintf(head01,'%s %-6d %s %-6d %s %s %s %s %10.1f\n',...
    'R:',2,':',i-1,':','*',':','*',R2(i));
       end
end

fclose(head01); % Closing the file 


toc