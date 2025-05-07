function [Fbest,Lbest,BestChart]=CPO(N,max_it,testsistem,Case)

    SearchAgents_no = N;
    Max_iter = max_it ;
    [lb,ub,dim]=ogranicenja(testsistem);

    Manis_pos=zeros(1,dim);
    Manis_score=inf;

    Ant_pos = zeros(1,dim);
    Ant_score=inf;

    fitness = zeros(SearchAgents_no, 1);

    %Initialize the positions of search agents
    Positions=initialization(SearchAgents_no,dim,ub,lb);
    Convergence_curve=zeros(1,Max_iter);
    t=0;
    % main loop
    while t<Max_iter
        for i=1:size(Positions,1)
            Flag4ub=Positions(i,:)>ub;
            Flag4lb=Positions(i,:)<lb;
            Positions(i,:)=(Positions(i,:) .* (~(Flag4ub+Flag4lb))) + ub.*Flag4ub + lb.*Flag4lb;
            % Calculate objective function for each search agent
            fitness(i)=Fobj(Positions(i,:),testsistem,Case);
            if fitness(i) <= Manis_score
                Manis_score=fitness(i); % Update Manis Pentadactyla
                Manis_pos = Positions(i,:);
            end
            if fitness(i) > Manis_score && fitness(i) < Ant_score
                Ant_score = fitness(i); % Update Ant
                Ant_pos = Positions(i,:);
            end
        end
            r1 = (rand()+rand())/2;
            r2 = rand();
            % Aroma concentration factor
            Cm = Aroma_concentration(Max_iter);
            % Rapid decrease factor
            C1 = (2-((t*2)/Max_iter));
            % Aroma trajectory factor
            a = Aroma_trajectory(SearchAgents_no,0.6);
            % Levy step length
            Levy_Step_length = Levy(SearchAgents_no) ;
            for i=1:SearchAgents_no
                % Energy correction factor
                lamda = 0.1*rand();
                VO2 = 0.2*rand();
                % Fatigue index factor
                Fatigue = log(((t*pi)/Max_iter)+1);
                % Energy consumption factor
                E = exp(-lamda*VO2*t*(1 + Fatigue));
                l = randi([1, Max_iter]);
                r3 = rand();

                % Energy fluctuation factor
                A1 = lamda*(2*E*rand()-E);
                % Luring behavior
                if Cm(l)>=0.2 && r3<=0.5
                    % Attraction and Capture Stage
                    D_ant =abs(a*Ant_pos-Manis_pos);
                    New_Ant_pos = Positions(i,:) + Ant_pos-A1*D_ant;
                    % Movement and Feeding Stage
                    D_manis = abs(C1*New_Ant_pos-Positions(i,:))-Levy_Step_length(i).*(1-t/Max_iter);
                    New_Manis_pos = Positions(i,:) +Manis_pos-A1*D_manis;
                    % Positions are updated
                    Positions(i,:) = (New_Manis_pos + New_Ant_pos)./2  ...
                        + ((sin(New_Ant_pos*exp(((t)/Max_iter))))...
                        ./((4*pi).*tan(New_Manis_pos*exp(((t*4*pi.^2)/(Max_iter))))))...
                        .*r1*r2*rand;
                % Predation behavior
                elseif Cm(l) <= 0.7 || r3 > 0.5
                    % Search and Localization Stage
                    if Cm(l)>=0 && Cm(l)<0.3
                        D_manis = abs(Levy_Step_length(i)*Manis_pos-Positions(i,:));
                        New_Manis_pos = sin((C1).*Positions(i,:)+A1.*abs(Manis_pos-Levy_Step_length(i).*D_manis));
                        Positions(i,:) = New_Manis_pos.*C1;
                    % Rapid Approach Stage
                    elseif Cm(l)>=0.3 && Cm(l)
                        D_manis = abs(a*Manis_pos-Positions(i,:));
                        New_Manis_pos = (Positions(i,:) - A1*abs(Manis_pos-exp(-a).*(rand.*pi)*D_manis));
                        Positions(i,:) = New_Manis_pos.*C1;
                    % Digging and Feeding Stage
                    elseif Cm(l)>=0.6
                        D_manis = abs(C1*Manis_pos-Positions(i,:));
                        New_Manis_pos = (Positions(i,:) + A1*abs(Manis_pos-D_manis));
                        Positions(i,:) = New_Manis_pos.*C1;
                    end
                end
            end
            t=t+1;
            Convergence_curve(t)=Manis_score;
    end    % end while
    Fbest = Manis_score;
    Lbest = Manis_pos;
    BestChart = Convergence_curve;
end

% This function initialize the first population of search agents
function Positions=initialization(SearchAgents_no,dim,ub,lb)
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
end

function Cm = Aroma_concentration(Max_iter)
    Q=100;
    % Preallocate arrays
    sigma_y = zeros(1, Max_iter);
    sigma_z = zeros(1, Max_iter);
    M = zeros(1, Max_iter);
    for t = 1:Max_iter
        r1=rand();
        H = 0.5 * r1;%Eq.(11)
        r2=rand();
        u = 2 + r2;%Eq.(10)
        sigma_y(t) = 50-((10*t)/Max_iter);%Eq.(12)
        sigma_z(t) = sin((pi*t)/(Max_iter))+40*exp(-t/Max_iter)-10*log((pi*t)/(Max_iter));%Eq.(13)
        M(t) = (Q/(pi*u*sigma_y(t)*sigma_z(t)))*exp(-(H^2)/(2*(sigma_z(t))^2));%Eq.(9)
        
    end
    HH=fi(M);
    Cm=rescale(HH);%Eq.(14)
end

function Bs = Aroma_trajectory(N, Dc)
    dt = 1 / N; 
    x0 = 0; y0 = 0; z0 = 0;
    
    % Initialize random number generator
    rng('default'); % Reset to default seed for reproducibility
    
    % Generate random steps (Brownian motion increments)
    dWx = sqrt(2 * Dc * dt) * randn(1, N);
    dWy = sqrt(2 * Dc * dt) * randn(1, N);
    dWz = sqrt(2 * Dc * dt) * randn(1, N);
    
    % Compute trajectory
    x = cumsum([x0, dWx]);
    y = cumsum([y0, dWy]);
    z = cumsum([z0, dWz]);
    
    % Seed based on current time for randomness
    rng('shuffle', 'twister'); 
    
    % Randomly sample a point from the trajectory
    randomIndex = randi(N);
    randomPoint = [x(randomIndex), y(randomIndex), z(randomIndex)];
    Bs = norm(randomPoint); % Euclidean distance from origin
end
function s1 = Levy(dim)
    beta = 1.5;
    sigma = (gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u = randn(1, dim) * sigma;
    v = randn(1, dim);
    step = u ./ abs(v).^(1/beta);
    s1=step;
end