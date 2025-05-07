function [] = ORPDexe(algorithm, testSystem, Case, numAgents, maxIter, numRuns)
    % ORPD Solver: Optimal Reactive Power Dispatch using Metaheuristics
    % Inputs:
    %   algorithm   - Optimization algorithm (e.g., 'GSA', 'PSO', 'GWO')
    %   testSystem  - Power system case (e.g., 'ts_ieee30', 'ts_ieee118')
    %   objective   - Objective function ('Ploss' or 'VD')
    %   numAgents   - Population size (default: 25)
    %   maxIter     - Max iterations (default: 50)
    %   numRuns     - Number of independent runs (default: 1)
    
    % Set Defaults & Input Validation
    if nargin < 6, numRuns = 1;       end
    if nargin < 5, maxIter = 50;     end
    if nargin < 4, numAgents = 25;    end
    if nargin < 3, Case = 'Ploss'; end
    if nargin < 2, testSystem = 'ts_ieee30'; end
    if nargin < 1, algorithm = 'PSOGSA'; end
    
    % Algorithm Metadata
    algorithmNames = containers.Map(...
        {'PSOGSA', 'GWO', 'BBO', 'BSA','ABC','TLBO','BSLO','FFA'}, ...
        {
         'Hybrid PSO-GSA', 'Grey Wolf', 'Biogeography-Based', ...
         'Backtracking Search','Artifecial Bee Colone','Teaching Learning Based Optimization'...
         'blood Sucking Leach optimizer','Fire Fly Algorithm'
         });
    
    % Display Configuration
    fprintf('\n[ORPD Configuration]');
    fprintf('\nAlgorithm: %s', algorithmNames(algorithm));
    fprintf('\nTest System: IEEE %s-bus', extractAfter(testSystem, 'ts_ieee'));
    fprintf('\nObjective: %s', getObjectiveName(Case));
    fprintf('\nAgents: %d | Iterations: %d | Runs: %d\n', numAgents, maxIter, numRuns);
    TPR=0;format long;
    while TPR<numRuns
        TPR=TPR+1;
        tic;
        methodMap = struct(...
        'PSOGSA', @PSOGSA, ...
        'GWO', @GWO, ...
        'FFA', @FFA, ...
        'ABC', @ABC, ...
        'TLBO', @TLBO, ...
        'GSASQP', @GSASQP, ...
        'FSA', @FSA,...
        'BSLO', @BSLO );

        if isfield(methodMap, algorithm)
            [Fbest, Lbest, BestChart] = methodMap.(algorithm)(numAgents, maxIter, testSystem, Case);
        else
            error('Unknown optimization method: %s', algorithm);
        end
        toc;
        runTime = toc; 

        bestObjectiveValuesPerRun(TPR) = Fbest;         % Optimal objective function value for this run
        bestControlVariablesPerRun(TPR, :) = Lbest;     % Optimal control variables for this run
        objectiveHistoryPerRun(TPR, :) = BestChart;     % Objective function convergence history for this run
        computationTimePerRun(TPR) = runTime;           % Execution time for this run
    end
        
        [~,bestObjectiveValuesIndex]=min(bestObjectiveValuesPerRun);
        contorl_vars=bestControlVariablesPerRun(bestObjectiveValuesIndex,:);
        objectiveHistory=objectiveHistoryPerRun(bestObjectiveValuesIndex,:);
        % Plot convergence curve for the best run
        plot(objectiveHistory, '--k', 'LineWidth', 2);
        title('\fontsize{11}\bf ORPD');
        xlabel('\fontsize{11}\bf Iteration'); 
        ylabel('\fontsize{11}\bf Fobj');
        grid on;

        % Dynamic legend based on method
        legend(sprintf('\\fontsize{11}\\bf %s', algorithm));

    


    %st(contorl_vars,testSystem,Case,numAgents,maxIter,numRuns,bestObjectiveValuesPerRun,runTime);
    % Add these lines to export variables to the workspace
     assignin('base', 'bestObjectiveValuesPerRun', bestObjectiveValuesPerRun);           % Optimal objective function values across runs
     assignin('base', 'bestControlVariablesPerRun', bestControlVariablesPerRun);         % Best control variables for each run
     assignin('base', 'objectiveHistoryPerRun', objectiveHistoryPerRun);           % Convergence curves for all runs
     assignin('base', 'computationTimePerRun', computationTimePerRun);             % Computation time for each run
     assignin('base', 'objectiveHistory', objectiveHistory);         % Best convergence curve (for plotting)

end
    
    %Helper Functions
    function name = getObjectiveName(obj)
        % Maps objective codes to full names
        switch obj
            case 'Ploss', name = 'Power Loss Minimization';
            case 'VD', name = 'Voltage Deviation Minimization';
        end
    end