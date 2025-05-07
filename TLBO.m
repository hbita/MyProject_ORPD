function [Fbest,Lbest,BestChart]=TLBO(numAgents,maxIter,testSystem,Case)
    [VarMin,VarMax,dim]=constraints(testSystem);
    VarSize = [1 dim] ;
    %% TLBO prams
    MaxIt = maxIter;        % Maximum Number of Iterations
    nPop = numAgents;
    %% Initialization 
    % Empty Structure for Individuals
    empty_individual.Position = [];
    empty_individual.Cost = [];

    % Initialize Population Array
    pop = repmat(empty_individual, nPop, 1);

    % Initialize Best Solution
    BestSol.Cost = inf;

    % Initialize Population Members
    for i=1:nPop
        pop(i).Position = unifrnd(VarMin, VarMax, VarSize);
        pop(i).Cost = obj_fun(pop(i).Position,testSystem,Case);
        
        if pop(i).Cost < BestSol.Cost
            BestSol = pop(i);
        end
    end
    BestCosts = zeros(MaxIt,1);
    %% main loop
    for it=1:MaxIt
        % Calculate Population Mean
        Mean = 0;
        for i=1:nPop
            Mean = Mean + pop(i).Position;
        end
        Mean = Mean/nPop;
        % Select Teacher
        Teacher = pop(1);
        for i=2:nPop
            if pop(i).Cost < Teacher.Cost
                Teacher = pop(i);
            end
        end
        % Teacher Phase
        for i=1:nPop
            % Create Empty Solution
            newsol = empty_individual;
            
            % Teaching Factor
            TF = randi([1 2]);
            
            % Teaching (moving towards teacher)
            newsol.Position = pop(i).Position ...
                + rand(VarSize).*(Teacher.Position - TF*Mean);
            
            % Clipping
            newsol.Position = max(newsol.Position, VarMin);
            newsol.Position = min(newsol.Position, VarMax);
            
            % Evaluation
            newsol.Cost = obj_fun(newsol.Position,testSystem,Case);
            
            % Comparision
            if newsol.Cost < pop(i).Cost
                pop(i) = newsol;
                if pop(i).Cost < BestSol.Cost
                    BestSol = pop(i);
                end
            end
        end
        % Learner Phase
        for i=1:nPop
                
            A = 1:nPop;
            A(i)=[];
            j = A(randi(nPop-1));
            
            Step = pop(i).Position - pop(j).Position;
            if pop(j).Cost < pop(i).Cost
                Step = -Step;
            end
            
            % Create Empty Solution
            newsol = empty_individual;
            
            % Teaching (moving towards teacher)
            newsol.Position = pop(i).Position + rand(VarSize).*Step;
            
            % Clipping
            newsol.Position = max(newsol.Position, VarMin);
            newsol.Position = min(newsol.Position, VarMax);
            
            % Evaluation
            newsol.Cost = obj_fun(newsol.Position,testSystem,Case);
            
            % Comparision
            if newsol.Cost<pop(i).Cost
                pop(i) = newsol;
                if pop(i).Cost < BestSol.Cost
                    BestSol = pop(i);
                end
            end
        end

        BestCosts(it) = BestSol.Cost;
    end        %main end
    % Store Record for Current Iteration
    Fbest=BestCosts(MaxIt);
    Lbest=BestSol.Position;
    BestChart=BestCosts;