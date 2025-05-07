% Gravitational Search Algorithm.
%testsistem='ts_ieee30';
%Case='Ploss';
function [Fbest,Lbest,BestChart]=GSASQP(N,max_it,testsistem,Case)
global testsistem Case

%V:   Velocity.
%a:   Acceleration.
%M:   Mass.  Ma=Mp=Mi=M;
%dim: Dimension of the test function.
%N:   Number of agents.
%X:   Position of agents. dim-by-N matrix.
%R:   Distance between agents in search space.
%[low-up]: Allowable range for search space.
%Rnorm:  Norm 
%Rpower: Power of R.

 ElitistCheck=1; Rpower=1;
 min_flag=1; % 1: minimization, 0: maximization
 Rnorm=2; 
 
%get allowable range and dimension of the test function.
[low,up,dim]=ogranicenja(testsistem); 

%random initialization for agents.
X=initialization(dim,N,up,low); 

%create the best so far chart and average fitnesses chart.
BestChart=[];
%MeanChart=[];

V=zeros(N,dim);

for iteration=1:max_it
%     iteration
    
    %Checking allowable range. 
    X=space_bound(X,up,low); 

    %Evaluation of agents. 
    fitness=evaluateF(X); 
    
    if min_flag==1
    [best best_X]=min(fitness); %minimization.
    else
    [best best_X]=max(fitness); %maximization.
    end        
    
    if iteration==1
       Fbest=best;Lbest=X(best_X,:);
    end
    if min_flag==1
      if best<Fbest  %minimization.
       Fbest=best;Lbest=X(best_X,:);
      end
    else 
      if best>Fbest  %maximization
       Fbest=best;Lbest=X(best_X,:);
      end
    end
%%-------------------------------------------------------------------------
RateLS=0.95;
%RateLS=1;
%if iteration>max_it*0.50
if rand>RateLS % We decide to perform a Local Search.
    options=optimset('MaxFunEvals',1e4,'Display','off',...
    'algorithm','sqp','UseParallel','never');
    [Lbest,Fbest] = fmincon(@Fobj_sqp,Lbest,[],[],[],[],low,up,[],options);
%    iteration
end
%end
%%-------------------------------------------------------------------------       
BestChart=[BestChart Fbest];
%MeanChart=[MeanChart mean(fitness)];

%Calculation of M. 
[M]=massCalculation(fitness,min_flag); 

%Calculation of Gravitational constant.
G=Gconstant(iteration,max_it); 

%Calculation of accelaration in gravitational field. 
a=Gfield(M,X,G,Rnorm,Rpower,ElitistCheck,iteration,max_it);

%Agent movement. 
[X,V]=move(X,a,V);

%--------------------------------------------------------------------------
%fprintf('GSA| Iter:%3d -->  Fbest: %9.5f\n',iteration,Fbest);
%--------------------------------------------------------------------------
end %iteration

%FUNCTIONS
%==========================================================================
%This function Evaluates the agents. 
function   fitness=evaluateF(X)
[N,dim]=size(X);
for i=1:N 
    %L is the location of agent number 'i'
    L=X(i,:); 
    %calculation of objective function for agent number 'i'
    fitness(i)=Fobj_sqp(L);
end
return

function [M]=massCalculation(fit,min_flag);
%%%%here, make your own function of 'mass calculation'
Fmax=max(fit); Fmin=min(fit); Fmean=mean(fit); 
[i N]=size(fit);
if Fmax==Fmin
   M=ones(N,1);
else  
   if min_flag==1 %for minimization
      best=Fmin;worst=Fmax; 
   else %for maximization
      best=Fmax;worst=Fmin; 
   end  
   M=(fit-worst)./(best-worst); 
end
M=M./sum(M); %eq. 16.
return

% This function calculates Gravitational constant. 
function G=Gconstant(iteration,max_it)
%%%here, make your own function of 'G'
  alfa=10;
  G0=100;
  G=G0*exp(-alfa*iteration/max_it);
  return
  
  %This function calculates the accelaration of each agent in gravitational field. eq.7-10,21.
function a=Gfield(M,X,G,Rnorm,Rpower,ElitistCheck,iteration,max_it);
[N,dim]=size(X);
 final_per=2; %In the last iteration, only 2 percent of agents apply force to the others.
%%%%total force calculation
 if ElitistCheck==1
     kbest=final_per+(1-iteration/max_it)*(100-final_per); 
     kbest=round(N*kbest/100);
 else
     kbest=N; %eq.9.
 end
    [Ms ds]=sort(M,'descend');
 for i=1:N
     E(i,:)=zeros(1,dim);
     for ii=1:kbest
         j=ds(ii);
         if j~=i
            R=norm(X(i,:)-X(j,:),Rnorm); %Euclidian distanse.
         for k=1:dim 
             E(i,k)=E(i,k)+rand*(M(j))*((X(j,k)-X(i,k))/(R^Rpower+eps));
              %note that Mp(i)/Mi(i)=1
         end
         end
     end
 end

%%acceleration
a=E.*G; %note that Mp(i)/Mi(i)=1
return

%This function updates the velocity and position of agents.
function [X,V]=move(X,a,V)
%movement.
[N,dim]=size(X);
V=rand(N,dim).*V+a; 
X=X+V; 
return

%This function checks the search space boundaries for agents.
function  X=space_bound(X,up,low)
[N,dim]=size(X);
for i=1:N 
%     %%Agents that go out of the search space, are reinitialized randomly .
    Tp=X(i,:)>up;Tm=X(i,:)<low;X(i,:)=(X(i,:).*(~(Tp+Tm)))+((rand(1,dim).*(up-low)+low).*(Tp+Tm));
%     %%Agents that go out of the search space, are returned to the boundaries.
%         Tp=X(i,:)>up;Tm=X(i,:)<low;X(i,:)=(X(i,:).*(~(Tp+Tm)))+up.*Tp+low.*Tm;

end
return

%This function initializes the position of the agents in the search space, randomly.
function [X]=initialization(dim,N,up,down)
if size(up,2)==1
    X=rand(N,dim).*(up-down)+down;
end
if size(up,2)>1
    for i=1:dim
    high=up(i);low=down(i);
    X(:,i)=rand(N,1).*(high-low)+low;
    end
end
return
