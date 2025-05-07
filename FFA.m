function [Fbest,Lbest,BestChart]=FFA(numAgents,maxIter,testSystem,Case)
    alpha=0.5;       % Randomness 0--1 (highly random)
    betamn=0.2;     % minimum value of beta
    gamma=1;         % Absorption coefficient
    para=[numAgents maxIter alpha betamn gamma];
    %get allowable range and dimension of the test function.
    [Lb,Ub,dim]=constraints(testSystem);
    % Initial random guess
    u0=Lb+(Ub-Lb).*rand(1,dim);
    n=para(1);
    MaxGeneration=para(2);
    alpha=para(3); 
    betamin=para(4); 
    gamma=para(5);
    % Total number of function evaluations
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
    for k=1:MaxGeneration           % start iterations

        % This line of reducing alpha is optional
        alpha=alpha_new(alpha,MaxGeneration);
        
        % Evaluate new solutions (for all n fireflies)
        for i=1:n
            zn(i)=obj_fun(ns(i,:),testSystem,Case);
            Lightn(i)=zn(i);
        end
        
        % Ranking fireflies by their light intensity/objectives
        [Lightn,Index]=sort(zn);
        ns_tmp=ns;
        for i=1:n
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
        
    end                     %end of iterations
    Fbest=fbest;
    Lbest=nbest;
    BestChart=fbestVector;
end


% The initial locations of n fireflies
function [ns,Lightn]=init_ffa(n,d,Lb,Ub,u0)
    % if there are bounds/limits,
  if length(Lb)>0,
     for i=1:n,
     ns(i,:)=Lb+(Ub-Lb).*rand(1,d);
     end
  else
     % generate solutions around the random guess
     for i=1:n,
     ns(i,:)=u0+randn(1,d);
     end
  end
  
  % initial value before function evaluations
  Lightn=ones(n,1)*10^100;
end


% Move all fireflies toward brighter ones
function [ns]=ffa_move(n,d,ns,Lightn,nso,Lighto,...
    ~,~,alpha,betamin,gamma,Lb,Ub)
    % Scaling of the system
    scale=abs(Ub-Lb);

    % Updating fireflies
    for i=1:n
        % The attractiveness parameter beta=exp(-gamma*r)
        for j=1:n
            r=sqrt(sum((ns(i,:)-ns(j,:)).^2));
            % Update moves
            if Lightn(i)>Lighto(j), % Brighter and more attractive
                 beta0=1; beta=(beta0-betamin)*exp(-gamma*r.^2)+betamin;
                 tmpf=alpha.*(rand(1,d)-0.5).*scale;
                 ns(i,:)=ns(i,:).*(1-beta)+nso(j,:).*beta+tmpf;
            end
        end % end for j

    end % end for i

    % Check if the updated solutions/locations are within limits
    [ns]=findlimits(n,ns,Lb,Ub);
    end

    % This function is optional, as it is not in the original FA
    % The idea to reduce randomness is to increase the convergence,
    % however, if you reduce randomness too quickly, then premature
    % convergence can occur. So use with care.
    function alpha=alpha_new(alpha,NGen)
        % alpha_n=alpha_0(1-delta)^NGen=10^(-4);
        % alpha_0=0.9
        delta=1-(10^(-4)/0.9)^(1/NGen);
        alpha=(1-delta)*alpha;
    end


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
end


% d-dimensional objective function
% function z=Fun(fhandle,nonhandle,u)
%      Objective
%      z=fhandle(u);
    
%      Apply nonlinear constraints by the penalty method
%      Z=f+sum_k=1^N lam_k g_k^2 *H(g_k) where lam_k >> 1
%      z=z+getnonlinear(nonhandle,u);
% end


% function Z=getnonlinear(nonhandle,u)
%     Z=0;
%     % Penalty constant >> 1
%     lam=10^15; lameq=10^15;
%     % Get nonlinear constraints
%     [g,geq]=nonhandle(u);
    
%     % Apply inequality constraints as a penalty function
%     for k=1:length(g)
%         Z=Z+ lam*g(k)^2*getH(g(k));
%     end
%     % Apply equality constraints (when geq=[], length->0)
%     for k=1:length(geq)
%        Z=Z+lameq*geq(k)^2*geteqH(geq(k));
%     end
% end

% Test if inequalities hold
% H(g) which is something like an index function
% function H=getH(g)
%     if g<=0
%         H=0;
%     else
%         H=1;
%     end
% end
    

% % Test if equalities hold
% function H=geteqH(g)
%     if g==0
%         H=0;
%     else
%         H=1;
%     end
% end