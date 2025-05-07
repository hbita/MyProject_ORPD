function [Fbest,Lbest,BestChart] = BSLO(N ,max_it,testsistem,Case)
%% problem definition
[lb,ub,dim]=ogranicenja(testsistem);
%%BSLO param
maxIter=max_it;
populationSize=N;
% Initialize population
leeches = rand(populationSize, dim) .* (ub - lb) + lb;
fitness = zeros(populationSize, 1);
% Evaluate initial fitness
for i = 1:populationSize
    fitness(i) = Fobj(leeches(i,:),testsistem,Case);
end

% Track the best solution
[bestFitness, bestIdx] = min(fitness);
bestSolution = leeches(bestIdx, :);

% Initialize convergence curve
convergenceCurve = zeros(1, maxIter);
convergenceCurve(1) = bestFitness;
% Main optimization loop
for iter = 2:maxIter
    for i = 1:populationSize
        % Update position (Exploitation: Move towards best solution)
        r = rand();
        if r < 0.7 % Exploitation probability
            leeches(i, :) = leeches(i, :) + r * (bestSolution - leeches(i, :));
        else % Exploration: Random movement
            leeches(i, :) = rand(1, dim) .* (ub - lb) + lb;
        end

        % Apply boundary constraints
        leeches(i, :) = max(min(leeches(i, :), ub), lb);

        % Evaluate fitness of the updated position
        newFitness = Fobj(leeches(i,:),testsistem,Case);

        % Update if better solution found
        if newFitness < fitness(i)
            fitness(i) = newFitness;

            % Update global best solution
            if newFitness < bestFitness
                bestFitness = newFitness;
                bestSolution = leeches(i, :);
            end
        end
    end
    % Update convergence curve
    convergenceCurve(iter) = bestFitness;

    % Display iteration info
    % fprintf('Iteration %d: Best Fitness = %.6f\n', iter, bestFitness);
end
Fbest=bestFitness;
Lbest=bestSolution;
BestChart=convergenceCurve;
end
