function [Y, branchAdmittance, shuntAdmittance] = Ybus(branch_data, bus_data,tap_ratio)
    % This function calculates the bus admittance matrix (Ybus) for a power system network
    
    [numBranches, ~] = size(branch_data);
    fromNode = branch_data(:,1);           % Starting bus of each branch
    toNode = branch_data(:,2);             % Ending bus of each branch
    branchType = branch_data(:,3);         % Type of each branch (1=line, 2=transformer)
    branchResistance = branch_data(:,4);   % Branch resistance (R)
    branchReactance = branch_data(:,5);    % Branch reactance (X)
    branchSusceptance = branch_data(:,6);  % Branch shunt susceptance (B)
    
    branchShuntAdmittance = 1i * branchSusceptance; % Shunt admittance of each branch
    
    % Branch impedance (R + jX)
    branchImpedance = branchResistance + 1i * branchReactance;
    
    [numNodes, ~] = size(bus_data);
    % Shunt admittance at each node (capacitors, reactors, etc.)
    nodeShuntAdmittance = bus_data(1:numNodes, 9);
    nodeShuntAdmittance = diag(1i * nodeShuntAdmittance); 
    
    % Initialize matrices
    impedanceMatrix = zeros(numNodes);
    branchAdmittance = zeros(numNodes);
    shuntAdmittance = zeros(numNodes);
    
    for k = 1:numBranches
        % Calculate elements of impedance and admittance matrices
        impedanceMatrix(fromNode(k), toNode(k)) = branchImpedance(k);
        impedanceMatrix(toNode(k), fromNode(k)) = branchImpedance(k);
        
        branchAdmittance(fromNode(k), toNode(k)) = 1/impedanceMatrix(fromNode(k), toNode(k))/tap_ratio(k);
        branchAdmittance(toNode(k), fromNode(k)) = 1/impedanceMatrix(toNode(k), fromNode(k))/tap_ratio(k);
        
        if branchType(k) == 1  % Transmission line
            shuntAdmittance(fromNode(k), toNode(k)) = branchShuntAdmittance(k)/2;
            shuntAdmittance(toNode(k), fromNode(k)) = branchShuntAdmittance(k)/2;
        elseif branchType(k) == 2  % Transformer
            shuntAdmittance(fromNode(k), toNode(k)) = (1-tap_ratio(k))/tap_ratio(k)^2 * (1/branchImpedance(k));
            shuntAdmittance(toNode(k), fromNode(k)) = (tap_ratio(k)-1)/tap_ratio(k) * (1/branchImpedance(k));
        end
    end
    
    % Add node shunt admittance to the shunt admittance matrix
    shuntAdmittance = shuntAdmittance + nodeShuntAdmittance;
    
    numBuses = max(max(fromNode, toNode));
    Y = zeros(numBuses);
    
    % Calculate elements of the Ybus matrix
    for ii = 1:numBuses
        for jj = 1:numBuses
            if ii == jj
                % Diagonal elements: sum of all admittances connected to the bus
                Y(ii,jj) = sum(branchAdmittance(ii,:)) + sum(shuntAdmittance(ii,:));
            else
                % Off-diagonal elements: negative of branch admittance
                Y(ii,jj) = -branchAdmittance(ii,jj);
            end
        end
    end
end