function [objfun] =obj_fun (control_vars,testSystem,Case)
% =============================================
% CALCULATION OF THE MULTI-OBJECTIVE FUNCTION IN OPTIMAL POWER FLOW.
% Stott's decoupled method is used for power flow calculations.
% objfun - objective (multi-objective) function
% contorl_vars - vector of control variables;
% =============================================
%Test system setup
switch testSystem
    case 'ts_ieee30'
        [branch_data,bus_data,generator_data,transformer_data,shunt_comp_data,Vpqmin,Vpqmax,~]=ts_ieee30;
    case 'ts_ieee118'
        [branch_data,bus_data,generator_data,transformer_data,shunt_comp_data,Vpqmin,Vpqmax,~]=ts_ieee118;
end
%DATA extraction
%bus data matrix
[num_buses,~] = size(bus_data);
bus = bus_data(:,1);                % bus numbers
type = bus_data(:,2);               % bus types
V0 = bus_data(:,3);                 % initial voltage
teta0 = bus_data(:,4);              % initial voltage angle
Pg = bus_data(:,5);                 % active power injected at bus (p.u)
Qg = bus_data(:,6);                 % reavtive power injected at bus (p.u)
Pd = bus_data(:,7);                 % active power demad (p.u)
Qd = bus_data(:,8);                 % reactive power demad (p.u)
%--------------------------------------------------------------------------
tap_ratio = branch_data(:,7);       % transformer tap ratio
%--------------------------------------------------------------------------
%generator data matrix
[num_gen,~] = size(generator_data);
gen_bus = generator_data(:,1);               % Bus numbers where generators are connected
%Vgmin = generator_data(:,2);                 % Minimum allowable voltage (p.u.)
%Vgmax = generator_data(:,3);                 % Maximum allowable voltage (p.u.)
Pgmin = generator_data(:,4);                 % Minimum active power output (p.u.)
Pgmax = generator_data(:,5);                 % Maximum active power output (p.u.)
Qgmin = generator_data(:,6);                 % Minimum reactive power output (p.u.)
Qgmax = generator_data(:,7);                 % Maximum reactive power output (p.u.)
%a = generator_data(:,8);                     % Constant cost term ($/hr)
%b = generator_data(:,9);                     % Linear cost coefficient ($/MWhr)
%c = generator_data(:,10);                    % Quadratic cost coefficient ($/MW²hr)
%d = generator_data(:,11);                    % Valve-point loading coefficient (non-convex cost)
%e = generator_data(:,12);                    % Valve-point loading offset (non-convex cost)
%--------------------------------------------------------------------------
%transformer data matrix
[num_tran,~] = size(transformer_data);                    % Number of controllable transformers
tran_loc = transformer_data(:,1);                       % Branch numbers where transformers are located
%tap_ratio_min = transformer_data(:,2);                  % Minimum tap ratio (e.g., 0.9)
%tap_ratio_max = transformer_data(:,3);                  % Maximum tap ratio (e.g., 1.1)
%--------------------------------------------------------------------------
%shunt compensator data matrix
[num_shunt,~] =size(shunt_comp_data);                   % Number of Shunt compensators
shunt_loc= shunt_comp_data(:,1);                        % Bus numbers where compensators are connected
%Qgcmin = shunt_comp_data(:,2);                          % Minimum reactive power output (p.u.)
%Qgcmax = shunt_comp_data(:,3);                          % Maximum reactive power output (p.u.)
%--------------------------------------------------------------------------
% Extract bus data using vectorized operations
slack = bus(type == 0);  % Slack bus(es) (type 0)
%PV_buses = bus(type == 1);      % PV nodes (type 1)
PQ_buses = bus(type == 2);      % PQ nodes (type 2)

% Get counts (if needed)
num_PQ = numel(PQ_buses);
%um_PV = numel(PV_buses);
%num_slack = numel(slack);
for k = 1:num_gen 
    V0(gen_bus(k)) = control_vars(k) ;
end
for k=(num_gen + 1):(num_gen + num_tran)
    tap_ratio(tran_loc(k- num_gen)) = control_vars(k);
end
for k=(num_gen + num_tran + 1):(num_gen + num_tran + num_shunt)
    Qg(shunt_loc(k - num_gen - num_tran)) = control_vars(k);
end
% Power injection at nodes
Pinj = Pg -Pd ;
Qinj = Qg -Qd ;
[Y,branchAdmittance,shuntAdmittance]=Ybus(branch_data,bus_data,tap_ratio);
G=real(Y); B=imag(Y);
G_branch = real(branchAdmittance);  % Branch conductance (g)
B_branch = imag(branchAdmittance);  % Branch susceptance (b)
%G_shunt = real(shuntAdmittance);   % Shunt conductance (g0)
B_shunt = imag(shuntAdmittance);    % Shunt susceptance (b0)
B1 = B - B_shunt ;
for k = 1:num_buses
    if type(k) == 0
        B1(:,k)=[];
        B1(k,:)=[];
    end
end
B2rr = B - B_shunt ;
for k = 1:num_buses
    if type(k) == 0 || type(k) == 1
        B2rr(:,k) = 0;
        B2rr(k,:) = 0;
    end
    if type(k) == 2 
        index(k) = k ;
    end
end
index=nonzeros(index);
for k=1:length(index)
    B2r(k,:)=B2rr(index(k),:);
end
for k=1:length(index)
    B2(:,k)=B2r(:,index(k));
end
iter=0;
itermax=100;
epsilon=0.00001;
maxDeltaP=epsilon+1;
V=V0;
teta=teta0;
while (iter < itermax) && (maxDeltaP > epsilon) 
    iter=iter+1;
    tetanew=teta;Vnew=V;
    for k = 1: num_buses
        sum_I=0;
        for j = 1: num_buses
            sum_I = sum_I + Vnew(j) * (G(k,j) * cos(tetanew(k) - tetanew(j))...
                                    + B(k,j) * sin(tetanew(k) - tetanew(j)));
        end
        Delta_P (k) = Pinj(k) -Vnew(k) * sum_I ;
    end
    %Angle Update
    Delta_P(slack) = [];       % Remove slack bus
    Vnew(slack) = [];          % Remove slack bus
    tetanew(slack) = [];       % Remove slack bus
    teta (slack) = [];
    Delta_teta = - B1 \ (Delta_P' ./ Vnew);  % Solve Δθ
    teta = tetanew + Delta_teta;                % Update angles
    %Slack Bus Reinsertion
    teta_temp = zeros(num_buses,1);
    teta_temp(slack) = teta0(slack);                              % Slack angle fixed
    teta_temp(1:slack-1) = teta(1:slack-1);                       % Buses before slack
    teta_temp(slack+1:num_buses) = teta(slack:num_buses-1);       % Buses after slack
    teta = teta_temp;
    %Voltage Magnitude Update     
    Vnew_temp = zeros(num_buses,1);
    Vnew_temp(slack) = V0(slack);                                 % Slack voltage fixed
    Vnew_temp(1:slack-1) = Vnew(1:slack-1);                       % Buses before slack
    Vnew_temp(slack+1:num_buses) = Vnew(slack:num_buses-1);       % Buses after slack
    Vnew = Vnew_temp;
    for k=1:num_buses
        sum_I = 0;
        for j=1:num_buses
            sum_I = sum_I + Vnew(j)*(G(k,j)*sin(teta(k)-teta(j)) ...
                                 - B(k,j)*cos(teta(k)-teta(j)));
        end
        Delta_Q(k) = Qinj(k) - Vnew(k) * sum_I;    
    end
    B2i=-inv(B2);
    br1 = 0;
    for k = 1:num_buses
        if type(k) == 2                % PQ bus
          br1 = br1 + 1;
          V_pq(br1) = Vnew(k);         % Voltage magnitudes at PQ buses
          Delta_Qpq(br1) = Delta_Q(k); % Reactive power mismatches at PQ buses
        end
    end
    for kk = 1 : num_PQ
        Delta_V(kk) = B2i(kk, :) * (Delta_Qpq' ./ V_pq');
        V_pqnew(kk) = V_pq(kk) + Delta_V(kk);
    end
    br = 0;
    for k = 1:num_buses
        if type(k) == 2  % PQ bus
         br = br + 1;
         V(k) = V_pqnew(br);  % Update PQ bus voltages
        end
    end
    maxDeltaP=max(Delta_P);
end

%Power Injection Calculation
for k = 1:num_buses
    SUM_I_P = 0;
    SUM_I_Q = 0;
    for l = 1:num_buses
        SUM_I_P = SUM_I_P + V(l) * (G(k,l)*cos(teta(k)-teta(l)) + B(k,l)*sin(teta(k)-teta(l)));
        SUM_I_Q = SUM_I_Q + V(l) * (G(k,l)*sin(teta(k)-teta(l)) - B(k,l)*cos(teta(k)-teta(l)));
    end
    Pinj(k) = V(k) * SUM_I_P;  % Active power injection at bus k
    Qinj(k) = V(k) * SUM_I_Q;  % Reactive power injection at bus k
    Sinj(k) = Pinj(k) + Qinj(k)*sqrt(-1);  % Complex power injection
end
Pg = Pinj + Pd;  % Total active generation = injection + demand
Qg = Qinj + Qd;  % Total reactive generation = injection + demand
%Branch Power Flow Calculation
for i = 1:num_buses
    for j = 1:num_buses
        % Active power flow from bus i to j
        PF(i,j) = (V(i)^2)*G_branch (i,j) - V(i)*V(j)*(G_branch (i,j)*cos(teta(i)-teta(j))...
         + B_branch(i,j)*sin(teta(i)-teta(j)));
        
        % Reactive power flow (with shunt susceptance bgr0)
        QPF(i,j) = (-V(i)^2)*(B_branch(i,j) + B_shunt(i,j)) - V(i)*V(j)*(G_branch (i,j)*sin(teta(i)-teta(j))...
         - B_branch(i,j)*cos(teta(i)-teta(j)));
        
        % Reactive power flow (without shunt susceptance)
        QF_noshunt(i,j) = (-V(i)^2)*B_branch(i,j) - V(i)*V(j)*(G_branch (i,j)*sin(teta(i)-teta(j))...
         - B_branch(i,j)*cos(teta(i)-teta(j)));
    end
end
%System Loss Calculation
tot_Ploss = sum(Pinj);           % Total active power loss
tot_Qloss = sum(sum(QF_noshunt));      % Total reactive power loss (excluding shunts)
%sum of voltage deviations at all PQ buses
FdeltaVpq = 0;
for k = 1 : num_PQ
    FdeltaVpq = FdeltaVpq + abs(V(PQ_buses(k)) - 1);  % Sum of |V - 1| for PQ buses
end
%--------------------------------------------------------------------------
% Penalty weights for constraint enforcement in optimal power flow
lambdaPg = 1e6;         % Penalty weight for active power generation (Pg) constraint violations
lambdaVpq = 1e6;        % Penalty weight for voltage magnitude constraints (0.95 < V < 1.05 pu)
lambdaQg = 1e4;         % Penalty weight for reactive power generation (Qg) constraint violations
for k=1:num_gen
    if gen_bus(k) == slack
        if Pg(slack) < Pgmin(k)
          Pgsllim = Pgmin(k);
        elseif Pg(slack) > Pgmax(k)
             Pgsllim = Pgmax(k);
        else
             Pgsllim = Pg(slack);
         end
    end
end
Pf_Pgsl = lambdaPg * (Pg(slack)-Pgsllim)^2;
for k = 1 : num_PQ
    if V(PQ_buses(k)) < Vpqmin
        Vpqlim(k) = Vpqmin;
    elseif V(PQ_buses(k)) > Vpqmax
        Vpqlim(k) =Vpqmax;
    else 
        Vpqlim(k)=V(PQ_buses(k));
    end
    deltaVpq(k)=(V(PQ_buses(k))-Vpqlim(k))^2;
end
Pf_Vpq=lambdaVpq*sum(deltaVpq);
for k=1:num_gen
    if Qg(gen_bus(k))<Qgmin(k)
        Qglim(k)=Qgmin(k);
    elseif Qg(gen_bus(k))>Qgmax(k)
        Qglim(k)=Qgmax(k);
    else
        Qglim(k)=Qg(gen_bus(k));
    end
    deltaQg(k)=(Qg(gen_bus(k))-Qglim(k))^2;
end
Pf_Qg=lambdaQg*(sum(deltaQg));
Pf_Qg=lambdaQg*(sum(deltaQg));
% Objective function
switch Case
    case 'Ploss'              %Power loss minimization
    objfun = tot_Ploss + Pf_Pgsl + Pf_Vpq + Pf_Qg;
    case 'VD'               %Voltage profile improvement
    objfun = FdeltaVpq+Pf_Pgsl+Pf_Vpq+Pf_Qg;
end
return