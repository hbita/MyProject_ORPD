function[branch_data,bus_data,generator_data,transformer_data,shunt_comp_data,Vpqmin,Vpqmax,Sbase]=ts_ieee30

%All parameters and quantities are specified in per unit (p.u.)
%========================================================================================================================================
Sbase = 100 ; %MVA, Base power of the system
%----------------------------------------------------------------------------------------------------------------------------------------
%network_matrix: matrix where the system configuration and element parameters are specified
%----------------------------------------------------------------------------------------------------------------------------------------
%from_bus - starting node of the branch
%to_bus - ending node of the branch
%element_type - type of branch: 1 - line; 2 - transformer
%r - active resistance of the branch (p.u.); x - reactance of the branch (p.u.); b - susceptance of the branch (p.u.);
%tap_ratio- transformer tap ratio in p.u. (regulating tap ratio)
%Smax - maximum apparent power of the branch (thermal limit, MVA);   element_number - element serial number
               
          %|bus i  ->  bus j  |    type  |    r      |    x        |     b       |   tap_ratio   |   Smax     | branch_number  |
branch_data =  [
            1             2           1       0.0192       0.0575       0.0528            1          1.3            1
            1             3           1       0.0452       0.1852       0.0408            1          1.3            2
            2             4           1        0.057       0.1737       0.0368            1         0.65            3
            3             4           1       0.0132       0.0379       0.0084            1          1.3            4
            2             5           1       0.0472       0.1983       0.0418            1          1.3            5
            2             6           1       0.0581       0.1763       0.0374            1         0.65            6
            4             6           1       0.0119       0.0414        0.009            1          0.9            7
            5             7           1        0.046        0.116       0.0204            1          0.7            8
            6             7           1       0.0267        0.082        0.017            1          1.3            9
            6             8           1        0.012        0.042        0.009            1         0.32           10
            6             9           2            0        0.208            0        1.078         0.65           11
            6            10           2            0        0.556            0        1.069         0.32           12
            9            11           1            0        0.208            0            1         0.65           13
            9            10           1            0         0.11            0            1         0.65           14
            4            12           2            0        0.256            0        1.032         0.65           15
            12           13           1            0         0.14            0            1         0.65           16
            12           14           1       0.1231       0.2559            0            1         0.32           17
            12           15           1       0.0662       0.1304            0            1         0.32           18
            12           16           1       0.0945       0.1987            0            1         0.32           19
            14           15           1        0.221       0.1997            0            1         0.16           20
            16           17           1       0.0824       0.1932            0            1         0.16           21
            15           18           1        0.107       0.2185            0            1         0.16           22
            18           19           1       0.0639       0.1292            0            1         0.16           23
            19           20           1        0.034        0.068            0            1         0.32           24
            10           2            1       0.0936        0.209            0            1         0.32           25
            10           17           1       0.0324       0.0845            0            1         0.32           26
            10           21           1       0.0348       0.0749            0            1         0.32           27
            10           22           1       0.0727       0.1499            0            1         0.32           28
            21           22           1       0.0116       0.0236            0            1         0.32           29
            15           23           1          0.1        0.202            0            1         0.16           30
            22           24           1        0.115        0.179            0            1         0.16           31
            23           24           1        0.132         0.27            0            1         0.16           32
            24           25           1       0.1885       0.3292            0            1         0.16           33
            25           26           1       0.2544         0.38            0            1         0.16           34
            25           27           1       0.1093       0.2087            0            1         0.16           35
            28           27           2            0        0.396            0        1.068         0.65           36
            27           29           1       0.2198       0.4153            0            1         0.16           37
            27           30           1       0.3202       0.6027            0            1         0.16           38
            29           30           1       0.2399       0.4533            0            1         0.16           39
            8            28           1       0.0636          0.2       0.0428            1         0.32           40
            6            28           1       0.0169       0.0599        0.013            1         0.32           41];
%--------------------------------------------------------------------------------------------------------------------  
%buses_matrix: matrix specifying the quantities at the system nodes (buses)
%--------------------------------------------------------------------------------------------------------------------
%bus - bus index
%type - bus type: 1 - PV bus (voltage-controlled); 2 - PQ bus (load bus)
%voltage_magnitude - voltage magnitude (p.u.); phase_angle - voltage phase angle (radians);
%P_gen - active power generated (p.u.); Q_gen - reactive power generated (p.u.);
%P_load - active power consumed (p.u.); Q_load - reactive power consumed (p.u.);
%Yshunt - shunt admittance (p.u.) connected to the bus (e.g., capacitors/reactors modeled via admittance, not power)
%Vbase - base voltage in kV


 % | bus | type |  V   | teta |  Pg    |  Qg   | Pd     | Qd  | Yshunt | Vbase | 
bus_data=[
     1      0    1.050     0    0.0000    0      0       0       0      132
     2      1    1.040     0    0.8000    0      0.217   0.127   0      132
     3      2    1.000     0    0.0000    0      0.024   0.012   0      132
     4      2    1.000     0    0.0000    0      0.076   0.016   0      132
     5      1    1.010     0    0.5000    0      0.942   0.190   0      132
     6      2    1.000     0    0.0000    0      0.000   0.000   0      132
     7      2    1.000     0    0.0000    0      0.228   0.109   0      132
     8      1    1.010     0    0.2000    0      0.300   0.300   0      132
     9      2    1.000     0    0.0000    0      0.000   0.000   0      1
    10      2    1.000     0    0.0000    0      0.058   0.020   0      33
    11      1    1.050     0    0.2000    0      0.000   0.000   0      11
    12      2    1.000     0    0.0000    0      0.112   0.075   0      33
    13      1    1.050     0    0.2000    0      0.000   0.000   0      11
    14      2    1.000     0    0.0000    0      0.062   0.016   0      33
    15      2    1.000     0    0.0000    0      0.082   0.025   0      33
    16      2    1.000     0    0.0000    0      0.035   0.018   0      33
    17      2    1.000     0    0.0000    0      0.090   0.058   0      33
    18      2    1.000     0    0.0000    0      0.032   0.009   0      33
    19      2    1.000     0    0.0000    0      0.095   0.034   0      33
    20      2    1.000     0    0.0000    0      0.022   0.007   0      33
    21      2    1.000     0    0.0000    0      0.175   0.112   0      33
    22      2    1.000     0    0.0000    0      0.000   0.000   0      33
    23      2    1.000     0    0.0000    0      0.032   0.016   0      33
    24      2    1.000     0    0.0000    0      0.087   0.067   0      33
    25      2    1.000     0    0.0000    0      0.000   0.000   0      33
    26      2    1.000     0    0.0000    0      0.035   0.023   0      33
    27      2    1.000     0    0.0000    0      0.000   0.000   0      33
    28      2    1.000     0    0.0000    0      0.000   0.000   0      132   
    29      2    1.000     0    0.0000    0      0.024   0.009   0      33
    30      2    1.000     0    0.0000    0      0.106   0.019   0      33]; 
    %SPECIFICATION OF LIMITS FOR SYSTEM QUANTITIES
%=========================================================================================
    
%Voltage limits for PQ buses:    
%--------------------------------------- 
 %Vpqmin=0.95; Vpqmax=1.10; 
 Vpqmin=0.95; Vpqmax=1.10;


 %Generators: voltage and power limits, and cost function coefficients:
 %----------------------------------------------------------------------------------------
 %gen_bus - index of the bus where the generator is connected
 %Vgmin, Vgmax - minimum and maximum generator voltage limits, respectively
 %Pgmin, Pgmax - minimum and maximum active power output limits, respectively
 %Qgmin, Qgmax - minimum and maximum reactive power output limits, respectively
 %a, b, c - coefficients of the generator cost function Cg=a+b*Pg+c*Pg^2;
 %d, e - coefficients accounting for valve-point effect in the cost function Cg=a+b*Pg+c*Pg^2+|d*sin(e*(Pgmin-Pg))|

%| gen_bus| Vgmin | Vgmax | Pgmin | Pgmax | Qgmin | Qgmax |  a  |  b  |  c  |  d  |  e  |
generator_data=[
    1      0.95    1.10   0.50    2.50    -0.20   2.00     0    200   37.5   18    3.7
    2      0.95    1.10   0.20    0.80    -0.20   1.00     0    175   175    16    3.8
    5      0.95    1.10   0.15    0.50    -0.15   0.80     0    100   625    14    4.0
    8      0.95    1.10   0.10    0.35    -0.15   0.60     0    325   83.4   12    4.5
    11     0.95    1.10   0.10    0.30    -0.10   0.50     0    300   250    13    4.2
    13     0.95    1.10   0.12    0.40    -0.15   0.60     0    300   250    13.5  4.1];

%Transformers: Tap ratio limits for transformers:
%----------------------------------------------------------------------------------------
%tran_loc - serial number of the transformer (element) within the network_matrix
%tmin, tmax - minimum and maximum tap ratio values of the transformer, respectively
 
 %|tran_loc|  tmin | tmax | bus_i| bus_j|
 transformer_data=[
     11     0.90     1.10    6      9
     12     0.90     1.10    6     10
     15     0.90     1.10    4     12
     36     0.90     1.10    28    27];
%Reactive power compensators: Limits for reactive power compensation
%----------------------------------------------------------------------------------------
%bus - index of the bus where the reactive power compensator is connected
%Qc_min, Qc_max - minimum and maximum reactive power output of the compensator, respectively

%| tran_loc | tap_ratio_min | tap_ratio_max|
 shunt_comp_data=[
    10          0               0.05
    12          0               0.05
    15          0               0.05
    17          0               0.05
    20          0               0.05
    21          0               0.05
    23          0               0.05
    24          0               0.05
    29          0               0.05];
end 