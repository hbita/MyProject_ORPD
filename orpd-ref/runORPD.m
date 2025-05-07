function [] = runORPD(method ,testsistem ,Case ,N ,max_it,BRP)

%------------------------------------------------------------------------------------------------------------------------------------------------------
% inputs:
% testsistem: 'ts_wscc9', 'ts_ieee14', 'ts_ieee30', 'ts_ieee39', 'ts_ieee57', 'ts_ieee118', .. etc. 
% Case - objective function:  'Ploss' - Power loss minimization
%                             'VD'  - Voltage profile improvement (minimzation volt. dev. in PQ nodes)
%                             
% N - Number of agents.
% max_it -  Maximum number of iterations (T).
% BRP - number of runs. 
%-------------------------------------------------------------------------------------------------------------------------------------------------------
%outputs:
% Fbest: Best result. 
% Lbest: Best solution. The location of Fbest in search space.
% BestChart: The best so far Chart over iterations. 
%BRP=1;

%--------------------------------------------------------------------------
%% default arguments
if nargin < 6
    BRP = 1;                                    
    if nargin < 5
        max_it = 100;                          
        if nargin < 4
            N = 25;                           
            if nargin < 3
                Case = 'Ploss';                 
                if nargin < 2
                    testsistem = 'ts_ieee30';  
                    if nargin < 1
                        method = 'PSOGSA';         
                    end
                end
            end
        end
    end
end
%--------------------------------------------------------------------------------------------------------------------------------------
switch method
      case 'GSA'
      fprintf('\n  Run ORPD by Gravitational Search Algorithm (GSA)...');
      case 'PSO'
      fprintf('\n  Run ORPD by Particle Swarm Optimization (PSO)...');
      case 'PSOGSA'
      fprintf('\n  Run ORPD by Hybrid Particle Swarm Optimization and Gravitational Search Algorithm (PSOGSA)...');
      case 'BBO'
      fprintf('\n  Run ORPD by Biogeography-based optimization (BBO)...');
      case 'GWO'
      fprintf('\n  Run ORPD by Grey Wolf Optimizer (GWO)...');
      case 'BSA'
      fprintf('\n  Run ORPD by Backtracking Search Optimization Algorithm (BSA)...');
      case 'FFA'
      fprintf('\n  Run ORPD by FireFly Algorithm (FFA)...');
      case 'WDO'
      fprintf('\n  Run ORPD by Wind Driven Optimization (WDO)...');
      case 'ABC'
      fprintf('\n  Run ORPD by Artificial Bee Colony (ABC)...');
      case 'CS'
      fprintf('\n  Run ORPD by Cuckoo Search Algorithm (CS)...');
      case 'JANA'
      fprintf('\n  Run ORPD by Novi (JANA)...');
      case 'MSA'
      fprintf('\n  Run ORPD by Moth Swarm Algorithm (MSA)...');
      case 'TLBO'
      fprintf('\n  Run ORPD by Teaching-Learning-Based Optimization (TLBO)...'); 
      case 'GSASQP'
      fprintf('\n  Run ORPD by Hybrid GSA-SQP Algorithm (GSASQP)...'); 
      case 'BSLO'
      fprintf('\n  Run ORPD by BSLO...'); 
 end

      fprintf('\n  Number of agents (population size): %4d',N);
      fprintf('\n  Maximum number of iterations: %4d',max_it);
      fprintf('\n  Number of Runs: %4d',BRP);
disp(' ');
switch testsistem
      case 'ts_ieee30'
      fprintf('\n  IEEE 30-bus Test System');
      case 'ts_ieee118'
      fprintf('\n  IEEE 118-bus Test System');
end

switch Case
      case 'Ploss'              
      fprintf('\n  Objective Function (OF): Minimization of Power Loss');  
      case 'VD'              
      fprintf('\n  Objective Function (OF): Minimization Voltage Deviation in PQ Nodes (Voltage Profile Improvement)');
end
fprintf('\n  --------------------------------------------------------\n')
%--------------------------------------------------------------------------------------------------------------------------------------

TPR=0;format long;
while TPR<BRP
TPR=TPR+1;tic;
switch method
      case 'GSA'
      [Fbest,Lbest,BestChart]=GSA(N,max_it,testsistem,Case);
      case 'PSO'
      [Fbest,Lbest,BestChart]=PSO(N,max_it,testsistem,Case);
      case 'PSOGSA'
      [Fbest,Lbest,BestChart]=PSOGSA(N,max_it,testsistem,Case);
      case 'BBO'
      [Fbest,Lbest,BestChart]=BBO(N,max_it,testsistem,Case);
      case 'GWO'
      [Fbest,Lbest,BestChart]=GWO(N,max_it,testsistem,Case);
      case 'BSA'
      [Fbest,Lbest,BestChart]=bsa(N,max_it,testsistem,Case);
      case 'FFA'
      [Fbest,Lbest,BestChart]=FFA(N,max_it,testsistem,Case);
      case 'WDO'
      [Fbest,Lbest,BestChart]=WDO(N,max_it,testsistem,Case);
      case 'ABC'
      [Fbest,Lbest,BestChart]=ABC(N,max_it,testsistem,Case);
      case 'CS'
      [Fbest,Lbest,BestChart]=CS(N,max_it,testsistem,Case);
      case 'JANA'
      [Fbest,Lbest,BestChart]=JANA(N,max_it,testsistem,Case);
      case 'MSA'
      [Fbest,Lbest,BestChart]=MSA(N,max_it,testsistem,Case);
      case 'TLBO'
      [Fbest,Lbest,BestChart]=TLBO(N,max_it,testsistem,Case);
      case 'GSASQP'
      [Fbest,Lbest,BestChart]=GSASQP(N,max_it,testsistem,Case);
      case 'FSA'
      [Fbest,Lbest,BestChart]=FSA(N,max_it,testsistem,Case);
       case 'BSLO'
      [Fbest,Lbest,BestChart]=BSLO(N,max_it,testsistem,Case);
      case 'FGO'
      [Fbest,Lbest,BestChart]=FGO(N,max_it,testsistem,Case);
end

%--------------------------------------------------------------------------------------------------------------------------------------
toc;
FBEST(TPR)=Fbest;                 % Vector containing the optimal values of the objective function for each execution of the calculation (TPR = 1 : BRP)
LBEST(TPR,:)=Lbest;               % Matrix whose rows are vectors of optimal values of control variables for each execution of the calculation (TPR = 1 : BRP) 
BESTCHART(TPR,:)=BestChart;      % Matrix whose rows are vectors with values of the objective function for each iteration of individual executions of the calculation
TOC(TPR)=toc;                     % Vector containing the duration of individual calculations (TPR = 1 : BRP)
end
[~,indFbest]=min(FBEST);      %minimalna (najbolja) vrednsot obj. funkcije i indeks prolaza (proracuna) u kome je ostvareno najbolje resenje 
L=LBEST(indFbest,:);              %optimalno resenje (vektor upravljackih promenljivih)
BstCh=BESTCHART(indFbest,:);  
%--------------------------------------------------------------------------
    
plot(BstCh,'--k','LineWidth',2);           %grafik Fobj za najbolji run
title('\fontsize{11}\bf ORPD');
xlabel('\fontsize{11}\bf Iteration');ylabel('\fontsize{11}\bf Fobj');
switch method
      case 'GSA'
        legend('\fontsize{11}\bf GSA');
      case 'PSO'
         legend('\fontsize{11}\bf PSO');
      case 'PSOGSA'
      legend('\fontsize{11}\bf PSOGSA');
         case 'BBO'
      legend('\fontsize{11}\bf BBO');
      case 'GWO'
         legend('\fontsize{11}\bf GWO');
      case 'BSA'
         legend('\fontsize{11}\bf BSA');
      case 'FFA'
         legend('\fontsize{11}\bf FFA');
      case 'WDO'
        legend('\fontsize{11}\bf WDO');
      case 'ABC'
         legend('\fontsize{11}\bf ABC');
      case 'CS'
         legend('\fontsize{11}\bf CS');
      case 'JANA'
        legend('\fontsize{11}\bf JANA');
      case 'MSA'
         legend('\fontsize{11}\bf MSA');
      case 'TLBO'
        legend('TLBO');
      case 'GSASQP'
        legend('\fontsize{11}\bf GSASQP');
      case 'FSA'
         legend('\fontsize{11}\bf FSA');
      case 'BSLO'
         legend('\fontsize{11}\bf BSLO');
      case 'FGO'
         legend('\fontsize{11}\bf FGO');
end
grid on;

st(L,testsistem,Case,N,max_it,BRP,FBEST,TOC);
% nr(L,testsistem,Case,N,max_it,BRP,FBEST,TOC)
% Add these lines to export variables to the workspace
% assignin('base', 'FBEST', FBEST);         % Optimal objective function values across runs
% assignin('base', 'LBEST', LBEST);         % Best control variables for each run
% assignin('base', 'BESTCHART', BESTCHART); % Convergence curves for all runs
% assignin('base', 'TOC', TOC);             % Computation time for each run
% assignin('base', 'BstCh', BstCh);         % Best convergence curve (for plotting)
return
