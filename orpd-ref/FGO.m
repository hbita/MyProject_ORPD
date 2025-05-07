function [Fbest,Lbest,BestChart]=FGO(N,max_it,testsistem,Case)  %N,max_it,testsistem,Case
    [lb,ub,dim]=ogranicenja(testsistem);
    %%%%-------------------Definitions--------------------------%%
    Lbest=zeros(1,dim); % A vector to include the best-so-far solution
    Fbest=inf; % A Scalar variable to include the best-so-far score
    BestChart=zeros(1,max_it);
    %%-------------------Controlling parameters--------------------------%%
    M=0.6;      % Determines the tradeoff percent between exploration and exploitation operators.
    Ep=0.7;     % Determines the probability of environmental effect on hyphal growth
    R=0.9;      % Determines the speed of convergence to the best-so-far solution
    %%---------------Initialization----------------------%%
    S=initialization(N,dim,ub,lb); % Initialize the S of crested porcupines
    t=0;        % Function evaluation counter
    fit = zeros(N, 1);
    %%---------------------Evaluation-----------------------%%
    for i=1:N
        % fit(i)=feval(fobj, S(i,:)');
        fit(i) = Fobj(S(i,:)',testsistem,Case);
    end
    % Update the best-so-far solution
    [Fbest,index]=min(fit);
    Lbest=S(index,:);
    % A new vector to store the best-so-far position for each hyphae
    Sp=S;

    while t<max_it 
        % ----------------------------------------------------------------------------- %%
        if t <= M/2                         % Compute the the nutrient allocation according to (12)
            nutrients = rand(N);      % Allocate more randomly to encourage fluctuation exploitation
        else
            nutrients = fit;                % Exploitation phase: allocate based on fitness
        end
        nutrients = nutrients / sum(nutrients)+2*rand; % Normalize nutrient allocation according to (13)
        if rand<rand    % Hyphal tip growth behavior 
            for i=1:N
                a=randi(N);
                b=randi(N);
                c=randi(N);
                while a==i || a==b || c==b || c==a ||c==i ||b==i
                    a=randi(N);
                    b=randi(N);
                    c=randi(N);
                end
                p=(fit(i)-min(fit))/(max(fit)-min(fit)+eps);      %Compute p_i according to (23)
                Er=M+(1-t/(max_it)).*(1-M);                      % Compute Er according to (24)
                if p<Er
                    F=(fit(i)/(sum(fit))).*rand*(1-t/(max_it))^(1-t/(max_it)); % Calculate F according to (5) and (6)
                    E=exp(F);  % Calculate E according to (4)
                    r1=rand(1,dim);
                    r2=rand;
                    U1=r1<r2; % Binary vector
                    S(i,:) = (U1).*S(i,:)+(1-U1).*(S(i,:)+E.*(S(a,:)-S(b,:))); % Generate the new hyphal growth according to (9)
                else
                    Ec = (rand(1, dim) - 0.5) .* rand.*(S(a,:)-S(b,:)); % Compute the additional exploratory step using (17)
                    if rand<rand      % Hypha growing in the opposite direction of nutrient-rich areas %
                        De2 = rand(1,dim).* (S(i,:) - Lbest).*(rand(1,dim)>rand); % Compute De2 according to (16) %%
                        S(i,:) = S(i,:) + De2 .* nutrients(i)+Ec*(rand>rand);
                    else % Growth direction toward nutrient-rich area
                        De = rand.* (S(a, :) - S(i, :)) + rand(1,dim).* ...
                                    ((rand>rand*2-1).*Lbest - S(i, :)).*(rand()>R);         % Compute De according to (11) %%
                        S(i,:) = S(i,:) + De .* nutrients(i)+Ec*(rand>Ep);                  % Compute the new growth for the ith hyphal using (14)
                    end
                end
                % Return the search agents that exceed the search space's bounds
                for j=1:size(S,2)
                    if  S(i,:)>ub
                        S(i,:)=lb+rand*(ub-lb);
                    elseif  S(i,j)<lb
                        S(i,:)=lb+rand*(ub-lb);
                    end
                    
                end
                % Calculate the fitness value of the newly generated solution
                nF=Fobj(S(i,:)',testsistem,Case);
                % update Global & Local best solution
                if  fit(i)<nF
                    S(i,:)=Sp(i,:);    % Update local best solution
                else
                    Sp(i,:)=S(i,:);
                    fit(i)=nF;
                    % update Global best solution
                    if  fit(i)<= Fbest
                        Lbest=S(i,:);    % Update global best solution
                        Fbest=fit(i);
                    end
                end
                t=t+1; % Move to the next generation
                if t>max_it
                    break
                end
                BestChart(t)=Fbest;
            end % End for i
        end % End If
        if  t>max_it
            break;
        end
    end  %end while
end 
% This function initialize the first population of search agents
function Positions=initialization(N,dim,ub,lb)

    Boundary_no= length(ub); % numnber of boundaries
    
    % If the boundaries of all variables are equal and user enter a signle
    % number for both ub and lb
    if Boundary_no==1
         Positions=rand(N,dim).*(ub-lb)+lb;
    end
    
    % If each variable has a different lb and ub
    if Boundary_no>1
        for i=1:dim
            ub_i=ub(i);
            lb_i=lb(i);
            Positions(:,i)=rand(N,1).*(ub_i-lb_i)+lb_i;      
        end
    end
end