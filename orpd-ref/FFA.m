% ======================================================== 
% Files of the Matlab programs included in the book:       
% Xin-She Yang, Nature-Inspired Metaheuristic Algorithms,  
% Second Edition, Luniver Press, (2010).   www.luniver.com 
% ========================================================    

% -------------------------------------------------------- 
% Firefly Algorithm for constrained optimization using     
% for the design of a spring (benchmark)                    
% by Xin-She Yang (Cambridge University) Copyright @2009   
% -------------------------------------------------------- 

%function bestsolution=fa_mincon(para)
function [Fbest,Lbest,BestChart]=FFA(N,max_it,testsistem,Case)
% ------------------------------------------------
% para = [n N_iteration alpha betamin gamma]
alpha=0.5;       % Randomness 0--1 (highly random)
betamn=0.2;     % minimum value of beta
gamma=1;         % Absorption coefficient
%para=[N max_it alpha betamn gamma];
para=[N max_it alpha betamn gamma];
% ------------------------------------------------
%format long;
% n=number of fireflies
% MaxGeneration=number of pseudo time steps
%Fbest=fbest;
%Lbest=nbest;
%BestChart=fbestVector;
%help fa_mincon.m
% This demo uses the Firefly Algorithm to solve the
% [Spring Design Problem as described by Cagnina et al.,
% Informatica, vol. 32, 319-326 (2008). ]

% Simple bounds/limits
%disp('Solve the simple spring design problem ...');
%disp('Solve OPF problem by FireFly Algorithm ...');
%Lb=[0.05 0.25 2.0];
%Ub=[2.0 1.3 15.0];

%get allowable range and dimension of the test function.
[Lb,Ub,dim]=ogranicenja(testsistem);

% Initial random guess
u0=Lb+(Ub-Lb).*rand(1,dim);

%[u,fval,NumEval]=ffa_mincon(@Fobj,@constraint,u0,Lb,Ub,para);
%[fval]=Fobj(u0,testsistem,Case);
% Display results
%bestsolution=u;
%bestojb=fval;
%total_number_of_function_evaluations=NumEval;

%%% Put your own cost/objective function here --------%%%
%% Cost or Objective function
 %function z=cost(x)
%z=(2+x(3))*x(1)^2*x(2);

% Constrained optimization using penalty methods
% by changing f to F=f+ \sum lam_j*g^2_j*H_j(g_j)
% where H(g)=0 if g<=0 (true), =1 if g is false

%%% Put your own constraints here --------------------%%%
%function [g,geq]=constraint(x)
% All nonlinear inequality constraints should be here
% If no inequality constraint at all, simple use g=[];
%g(1)=1-x(2)^3*x(3)/(71785*x(1)^4);
% There was a typo in Cagnina et al.'s paper, 
% the factor should 71785 insteady of 7178 !     
%tmpf=(4*x(2)^2-x(1)*x(2))/(12566*(x(2)*x(1)^3-x(1)^4));
%g(2)=tmpf+1/(5108*x(1)^2)-1;
%g(3)=1-140.45*x(1)/(x(2)^2*x(3));
%g(4)=x(1)+x(2)-1.5;
%g=[];
% all nonlinear equality constraints should be here
% If no equality constraint at all, put geq=[] as follows
%geq=[];

%%% End of the part to be modified -------------------%%%

%%% --------------------------------------------------%%%
%%% Do not modify the following codes unless you want %%%
%%% to improve its performance etc                    %%%
% -------------------------------------------------------
% ===Start of the Firefly Algorithm Implementation ======
% Inputs: fhandle => @cost (your own cost function,
%                   can be an external file  )
%     nonhandle => @constraint, all nonlinear constraints
%                   can be an external file or a function
%         Lb = lower bounds/limits
%         Ub = upper bounds/limits
%   para == optional (to control the Firefly algorithm)
% Outputs: nbest   = the best solution found so far
%          fbest   = the best objective value
%      NumEval = number of evaluations: n*MaxGeneration
% Optional:
% The alpha can be reduced (as to reduce the randomness)
% ---------------------------------------------------------
%fhandle=@Fobj;
% Start FA
%function [nbest,fbest,NumEval]...
%           =ffa_mincon(fhandle,nonhandle,u0, Lb, Ub, para)
% Check input parameters (otherwise set as default values)
%if nargin<6, para=[20 50 0.25 0.20 1]; end
%if nargin<5, Ub=[]; end
%if nargin<4, Lb=[]; end
%if nargin<3,
%disp('Usuage: FA_mincon(@cost, @constraint,u0,Lb,Ub,para)');
%end

% n=number of fireflies
% MaxGeneration=number of pseudo time steps
% ------------------------------------------------
% alpha=0.25;      % Randomness 0--1 (highly random)
% betamn=0.20;     % minimum value of beta
% gamma=1;         % Absorption coefficient
% ------------------------------------------------
n=para(1);  MaxGeneration=para(2);
alpha=para(3); betamin=para(4); gamma=para(5);

% Total number of function evaluations
NumEval=n*MaxGeneration;

% Check if the upper bound & lower bound are the same size
if length(Lb) ~=length(Ub)
    disp('Simple bounds/limits are improper!');
    return
end

% Calcualte dimension
d=dim;

% Initial values of an array
zn=ones(n,1)*10^100;
% ------------------------------------------------
% generating the initial locations of n fireflies
[ns,Lightn]=init_ffa(n,d,Lb,Ub,u0);

% Iterations or pseudo time marching
for k=1:MaxGeneration,     %%%%% start iterations

% This line of reducing alpha is optional
 alpha=alpha_new(alpha,MaxGeneration);

% Evaluate new solutions (for all n fireflies)
for i=1:n,
   %zn(i)=Fun(fhandle,nonhandle,ns(i,:));
   zn(i)=Fobj(ns(i,:),testsistem,Case);
   Lightn(i)=zn(i);
end

% Ranking fireflies by their light intensity/objectives
[Lightn,Index]=sort(zn);
ns_tmp=ns;
for i=1:n,
 ns(i,:)=ns_tmp(Index(i),:);
end

%% Find the current best
nso=ns; Lighto=Lightn;
nbest=ns(1,:); Lightbest=Lightn(1);

% For output only
fbest=Lightbest;

% Move all fireflies to the better locations
[ns]=ffa_move(n,d,ns,Lightn,nso,Lighto,nbest,...
      Lightbest,alpha,betamin,gamma,Lb,Ub);

fbestVector(k)=fbest;

%--------------------------------------------------------------------------
%fprintf('FFA|%5.0f -----> %9.5f\n',k,fbest);
%--------------------------------------------------------------------------

end   %%%%% end of iterations
Fbest=fbest;
Lbest=nbest;
BestChart=fbestVector;
%% initialize FFA plot
%plot(fbestVector,'--k','linewidth',2)
%title(['\fontsize{11}\bf OPF']);
%xlabel('\fontsize{11}\bf Iteration');ylabel('\fontsize{11}\bf Fobj');
%legend('\fontsize{11}\bf FFA',1);grid on;
% -------------------------------------------------------
% ----- All the subfunctions are listed here ------------
% The initial locations of n fireflies
function [ns,Lightn]=init_ffa(n,d,Lb,Ub,u0)
  % if there are bounds/limits,
if ~isempty(Lb)
   for i=1:n
   ns(i,:)=Lb+(Ub-Lb).*rand(1,d);
   end
else
   % generate solutions around the random guess
   for i=1:n
   ns(i,:)=u0+randn(1,d);
   end
end

% initial value before function evaluations
Lightn=ones(n,1)*10^100;

% Move all fireflies toward brighter ones
function [ns]=ffa_move(n,d,ns,Lightn,nso,Lighto,...
             nbest,Lightbest,alpha,betamin,gamma,Lb,Ub)
% Scaling of the system
scale=abs(Ub-Lb);

% Updating fireflies
for i=1:n
% The attractiveness parameter beta=exp(-gamma*r)
   for j=1:n
      r=sqrt(sum((ns(i,:)-ns(j,:)).^2));
      % Update moves
        if Lightn(i)>Lighto(j) % Brighter and more attractive
           beta0=1; beta=(beta0-betamin)*exp(-gamma*r.^2)+betamin;
           tmpf=alpha.*(rand(1,d)-0.5).*scale;
           ns(i,:)=ns(i,:).*(1-beta)+nso(j,:).*beta+tmpf;
        end
   end % end for j

end % end for i

% Check if the updated solutions/locations are within limits
[ns]=findlimits(n,ns,Lb,Ub);

% This function is optional, as it is not in the original FA
% The idea to reduce randomness is to increase the convergence,
% however, if you reduce randomness too quickly, then premature
% convergence can occur. So use with care.
function alpha=alpha_new(alpha,NGen)
% alpha_n=alpha_0(1-delta)^NGen=10^(-4);
% alpha_0=0.9
delta=1-(10^(-4)/0.9)^(1/NGen);
alpha=(1-delta)*alpha;

% Make sure the fireflies are within the bounds/limits
function [ns]=findlimits(n,ns,Lb,Ub)
for i=1:n,
     % Apply the lower bound
  ns_tmp=ns(i,:);
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I);

  % Apply the upper bounds
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
  % Update this new move
  ns(i,:)=ns_tmp;
end

% -----------------------------------------
% % d-dimensional objective function
function z=Fun(fhandle,nonhandle,u)
% Objective
z=fhandle(u);

% Apply nonlinear constraints by the penalty method
% Z=f+sum_k=1^N lam_k g_k^2 *H(g_k) where lam_k >> 1
z=z+getnonlinear(nonhandle,u);

function Z=getnonlinear(nonhandle,u)
Z=0;
% Penalty constant >> 1
lam=10^15; lameq=10^15;
% Get nonlinear constraints
[g,geq]=nonhandle(u);

% Apply inequality constraints as a penalty function
for k=1:length(g)
    Z=Z+ lam*g(k)^2*getH(g(k));
end
% Apply equality constraints (when geq=[], length->0)
for k=1:length(geq)
   Z=Z+lameq*geq(k)^2*geteqH(geq(k));
end

% Test if inequalities hold
% H(g) which is something like an index function
function H=getH(g)
if g<=0,
    H=0;
else
    H=1;
end

% Test if equalities hold
function H=geteqH(g)
if g==0,
    H=0;
else
    H=1;
end
%% ==== End of Firefly Algorithm implementation ======

