%  Grey Wolf Optimizer (GWO) source codes version 1.1               %
%                                                                   %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%                                                                   %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%       Homepage: http://www.alimirjalili.com/GWO.html              %
%                                                                   %
%   Main paper: S. Mirjalili, S. M. Mirjalili, A. Lewis             %
%               Grey Wolf Optimizer, Advances in Engineering        %
%               Software, Volume 69, March 2014, Pages 46-61,       %
%               http://dx.doi.org/10.1016/j.advengsoft.2013.12.007  %
%                                                                   %
% Grey Wolf Optimizer
%function [Alpha_score,Alpha_pos,Convergence_curve]=GWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj,handles,Value)
function [Fbest,Lbest,BestChart]=GWO(N,max_it,testsistem,Case)

%get allowable range and dimension of the test function.
[lb,ub,dim]=ogranicenja(testsistem); 

% initialize alpha, beta, and delta_pos
Lbest=zeros(1,dim);
Fbest=inf; %change this to -inf for maximization problems

Beta_pos=zeros(1,dim);
Beta_score=inf; %change this to -inf for maximization problems

Delta_pos=zeros(1,dim);
Delta_score=inf; %change this to -inf for maximization problems

%Initialize the positions of search agents
Positions=initializationGWO(N,dim,ub,lb);

%Convergence_curve=zeros(1,Max_iter);

l=0;% Loop counter

% Main loop
while l<max_it
    for i=1:size(Positions,1)             
        
        % Calculate objective function for each search agent
        fitness=Fobj(Positions(i,:),testsistem,Case);
        All_fitness(1,i)=fitness;
        
        % Update Alpha, Beta, and Delta
        if fitness<Fbest 
            Fbest=fitness; % Update alpha
            Lbest=Positions(i,:);
        end
        
        if fitness>Fbest && fitness<Beta_score 
            Beta_score=fitness; % Update beta
            Beta_pos=Positions(i,:);
        end
        
        if fitness>Fbest && fitness>Beta_score && fitness<Delta_score 
            Delta_score=fitness; % Update delta
            Delta_pos=Positions(i,:);
        end
    end
    
    
    a=2-l*((2)/max_it); % a decreases linearly fron 2 to 0
    
    % Update the Position of search agents including omegas
    for i=1:size(Positions,1)
        for j=1:size(Positions,2)     
                       
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand(); % r2 is a random number in [0,1]
            
            A1=2*a*r1-a; % Equation (3.3)
            C1=2*r2; % Equation (3.4)
            
            D_alpha=abs(C1*Lbest(j)-Positions(i,j)); % Equation (3.5)-part 1
            X1=Lbest(j)-A1*D_alpha; % Equation (3.6)-part 1
                       
            r1=rand();
            r2=rand();
            
            A2=2*a*r1-a; % Equation (3.3)
            C2=2*r2; % Equation (3.4)
            
            D_beta=abs(C2*Beta_pos(j)-Positions(i,j)); % Equation (3.5)-part 2
            X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2       
            
            r1=rand();
            r2=rand(); 
            
            A3=2*a*r1-a; % Equation (3.3)
            C3=2*r2; % Equation (3.4)
            
            D_delta=abs(C3*Delta_pos(j)-Positions(i,j)); % Equation (3.5)-part 3
            X3=Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3             
            
            Positions(i,j)=(X1+X2+X3)/3;% Equation (3.7)
            
        end
        
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
    end
%--------------------------------------------------------------------------
%fprintf('GWO|%5.0f -----> %9.5f\n',l,Fbest);
%--------------------------------------------------------------------------
    l=l+1;    
    BestChart(l)=Fbest;
    
    
    %if l>1
    %    line([l-1 l], [Convergence_curve(l-1) Convergence_curve(l)],'Color','k')
    %    xlabel('Iteration');
    %    ylabel('Best score obtained so far');        
    %    drawnow
    % end
 
    
    %set(handles.itertext,'String', ['The current iteration is ', num2str(l)])
    %set(handles.optimumtext,'String', ['The current optimal value is ', num2str(Alpha_score)])
    %if Value==1
    %    hold on
    %    scatter(l*ones(1,SearchAgents_no),All_fitness,'.','k')
    %end

end
%plot(BestChart,'--k','LineWidth',2)
%title(['\fontsize{11}\bf OPF']);
%xlabel('\fontsize{11}\bf Iteration');ylabel('\fontsize{11}\bf Fobj');
%legend('\fontsize{11}\bf GWO',1);grid on;

function Positions=initializationGWO(SearchAgents_no,dim,ub,lb)

Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    Positions=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end
return
