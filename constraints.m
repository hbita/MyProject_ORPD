function [down,up,dim]=constraints(testSystem)
    %If lower bounds of dimensions are the same, then 'down' is a value.
    %Otherwise, 'down' is a vector that shows the lower bound of each dimension.
    %This is also true for upper bounds of dimensions
    switch testSystem
        case 'ts_ieee30'
            [~,bus_data,generator_data,transformer_data,shunt_comp_data,~,~,~]=ts_ieee30;
        case 'ts_ieee118'
            [~,bus_data,generator_data,transformer_data,shunt_comp_data,~,~,~]=ts_ieee118;
    end
    %Extract Data from Loaded Variables
    [num_buses,~] = size(bus_data);
    bus = bus_data(:,1);                                    % bus numbers
    type = bus_data(:,2);    
    Vgmin = generator_data(:,2);                            % Minimum allowable voltage (p.u.)
    Vgmax = generator_data(:,3);   
    tap_ratio_min = transformer_data(:,2);                  % Minimum tap ratio (e.g., 0.9)
    tap_ratio_max = transformer_data(:,3);  
    Qgcmin = shunt_comp_data(:,2);                          % Minimum reactive power output (p.u.)
    Qgcmax = shunt_comp_data(:,3);
    %Identify Slack, PV, and PQ Nodes
    % Precompute the number of PQ and PV nodes
    num_pq = sum(type == 2);  % Count PQ nodes (type 2)
    num_pv = sum(type == 1);  % Count PV nodes (type 1)
    PQ = zeros(1, num_pq);  
    PV = zeros(1, num_pv); 
    %Reset counters
    num_PQ=0;
    num_PV=0;
    for k=1:num_buses
        if type(k) == 0
            slack=bus(k);
        elseif type(k) == 1
            num_PV = num_PV + 1;
            PV(num_PV) = bus(k);
        elseif type(k) == 2
            num_PQ = num_PQ +1 ;
            PQ(num_PQ) = bus(k) ;
        end
    end
    %Pgming = Pgmin;  % Copy of minimum active power limits
    %Pgmaxg = Pgmax;  % Copy of maximum active power limits
    %for k=1:num_buses
    %    if type(k) == slack
    %        Pgming(k) = [];  % Remove slack bus Pmin entry
    %        Pgmaxg(k) = [];  % Remove slack bus Pmax entry
    %    end
    %end

    %Construct Constraint Vectors
    down=[Vgmin;tap_ratio_min;Qgcmin]' ; %(lower bounds) 
    up=[Vgmax;tap_ratio_max;Qgcmax]' ;   %(upper bounds)
    dim=length(down);
end
