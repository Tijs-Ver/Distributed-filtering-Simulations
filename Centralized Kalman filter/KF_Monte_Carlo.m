clc; clear; close all; 


%% Initilization
sampletime = 5;                                                             %% Sample time [s]
Simulation_Time = 100;                                                      %% Simulation duration [s]
Time_Vector = 0:sampletime:Simulation_Time;                                 %% Create time vector with all values
Simulation_Steps = length(Time_Vector);                                     %% Number of time steps
Forcing_known_mode = false;
N_sim = 500;                                                                %% Number of simulations



%% Element properties                                                       %% Example material copper
rod_density = 8960;                                                         %% Unit [kg/m^3]
Thermal_diffusivity = 1.11e-4;                                              %% Unit [m^2/s]
specific_heat = 385;                                                        %% Unit [J/(kg * K)]

%% Geometric properties
L = 1;                                                                      %% Rod length [m]
A_cross = 1e-3;                                                             %% Cross-sectional area rod (Note this should be a lot smaller then the rod length, since only 1D conduction is considered)




%% Ground truth

%% Discretisation
N_Element_fine = 1000;                                                      %% Number of element in ground truth (Note: This value should be significantly larger than "N_Element_model" later on)
N_fine = N_Element_fine + 1;                                                %% Number of nodes approximating the true solution                                                          
N_Element_model = 100;                                                      %% Number of elements in model
N_model = N_Element_model + 1;                                              %% Number of nodes in model

Models_ratio = N_Element_fine/N_Element_model;                              %% Ratio between the ground truth approximation and model

%% Geometric properties ground truth
x = linspace(0, L, N_fine);                                                 %% Spatial positions[m]
dx = L / (N_fine - 1);                                                      %% Distance between nodes/Element length [m]
Volume_element_fine = dx * A_cross;                                         %% Unit [m^3]
Volume_element_model = Volume_element_fine * Models_ratio;                  %% Unit [m^3]

%% Gemetric properties model
x_model = linspace(0, L, N_model);                                          %% Spatial position along the model
dx_model = L/(N_model-1);                                                   %% Model element length


%% Initial condition and boundary condition
T_0 = ones(N_fine, 1) * 300;                                                %% Initial temperature rod [K]
Boundary_temperature= [300; 300];                                           %% Temperature the ends of the rods are held at [K] [left end; right end]

%Computing stiffness and mass matrix
Element_mass_matrix = dx/6 * [2,1;
                             1, 2];

Element_stiffness_matrix = Thermal_diffusivity/dx * [1, -1;
                                                    -1, 1];

Mass_matrix = zeros(N_fine,N_fine);

Stiffness_matrix = zeros(N_fine,N_fine);

for e = 1:N_fine-1                                                          %% Sum over all elements

    n1 = e;                                                                 %% Left node of the element
    n2 = e+1;                                                               %% Right node of the element

    %% Assembly
    Mass_matrix([n1 n2],[n1 n2]) = Mass_matrix([n1 n2],[n1 n2]) + Element_mass_matrix;
    Stiffness_matrix([n1 n2],[n1 n2]) = Stiffness_matrix([n1 n2],[n1 n2]) + Element_stiffness_matrix;
end



%% Forcing vector construction

loads = [
    0.30, -10;
    0.50, 30;
    0.70, -10];                                                             %% Loads defined as follows (Position along the rod [%], Heating/cooling [W])


forcing_fine = zeros(N_fine, 1);                                            %% Initialise forcing vector
forcing_model = zeros(N_model, 1);
for i = 1:size(loads, 1)                                                    %% loop over all the loads
    pos_ratio = loads(i, 1);                                                %% relative position along the rod
    power = loads(i, 2);                                                    %% Power at that element [W]
    
    Q = power / (rod_density * Volume_element_model * specific_heat);       %% Heating rate of the entire heated element [K/s]
    
    
    Element_forcing_vector = Q * (dx / 2) * [1; 
                                             1];                            %% [m*K/s]
    Element_forcing_vector_model = (Q * (dx_model / 2) * [1; 
                                             1]); 
    
    model_elem_idx = round(pos_ratio * N_Element_model);                    %% Corresponding element in the model (Elements 30, 50 and 70 in the example
    fine_start_elem = (model_elem_idx - 1) * round(Models_ratio) + 1;       %% First fine element part of the model element
    fine_end_elem   = model_elem_idx * round(Models_ratio);                 %% Last fine element part of the model element
    forcing_model([model_elem_idx, model_elem_idx+1]) = forcing_model([model_elem_idx, model_elem_idx+1]) + Element_forcing_vector_model;
    for e = fine_start_elem:fine_end_elem
        n1 = e;                                                             %% Left node of element
        n2 = e + 1;                                                         %% Right node of element
        
        %% Assembly
        forcing_fine([n1 n2]) = forcing_fine([n1 n2]) + Element_forcing_vector;
    end
end


%% Partioning


idx_b = [1, N_fine];                                                        %% Boundary indices
idx_i = 2:N_fine-1;                                                         %% Internal indices


M_ii = Mass_matrix(idx_i, idx_i);                                           %% M_{int,int} the effect of the temperature change of the internal points on the internal points                                           
K_ii = Stiffness_matrix(idx_i, idx_i);                                      %% K_{int,int} the effect of the temperature of the internal points on the internal points
K_ib = Stiffness_matrix(idx_i, idx_b);                                      %% K_{int,b} the effect of the temperature of the boundary on the internal points

f_i = forcing_fine(idx_i);                                                  %% F_{int} Forcing on the internal points

%% ODE15 Solving

RHS = f_i - K_ib * Boundary_temperature;                                    %% Right hand side of eq1.15 precalculated to save time

dT_dt = @(t, T_int) M_ii \ (RHS - K_ii * T_int);                            %% ODE 1.15 rewritten to isolate dT/dt


T0_int = T_0(idx_i);                                                        %% Initial temperature of internal nodes

[t_out, T_int_out] = ode15s(dT_dt, Time_Vector, T0_int);                    %% Call ODE15 function

True_temperature = [ones(size(t_out))*Boundary_temperature(1), T_int_out, ones(size(t_out))*Boundary_temperature(2)]'; %% Add boundary values

%% Interpolating ground truth to smaller grid

True_temperature_at_Model_Nodes = zeros(N_model, Simulation_Steps);         %% Preallocate sized down ground truth

for k = 1:Simulation_Steps                                                  %% Interpolate
    True_temperature_at_Model_Nodes(:, k) = interp1(x, True_temperature(:, k), x_model, 'linear')';
end


%% Model

%% Initial temperature and boundary conditions model

T_0 = ones(N_model, 1) * 300;                                               %% Initial temperature rod [K]
Boundary_temperature_model= [300; 300];                                     %% Temperature the ends of the rods are held at [K] [left end; right end]

%Computing stiffness and mass matrix
Element_mass_matrix_model = dx_model/6 * [2,1;
                             1, 2];

Element_stiffness_matrix_model = Thermal_diffusivity/dx_model * [1, -1;
                                                    -1, 1];

Mass_matrix_model = zeros(N_model,N_model);

Stiffness_matrix_model = zeros(N_model,N_model);

for e = 1:N_model-1                                                         %% Sum over all elements
    
    
    n1 = e;                                                                 %% Left node element
    n2 = e+1;                                                               %% Right node element

    %% Assembly
    Mass_matrix_model([n1 n2],[n1 n2]) = Mass_matrix_model([n1 n2],[n1 n2]) + Element_mass_matrix_model;
    Stiffness_matrix_model([n1 n2],[n1 n2]) = Stiffness_matrix_model([n1 n2],[n1 n2]) + Element_stiffness_matrix_model;
end


%% Partioning


idx_b_model = [1, N_model];                                                 %% Boundary indices
idx_i_model = 2:N_model-1;                                                  %% Internal indices


M_ii_model = Mass_matrix_model(idx_i_model, idx_i_model);                   %% M_{int,int} the effect of the temperature change of the internal points on the internal points
K_ii_model = Stiffness_matrix_model(idx_i_model, idx_i_model);              %% K_{int,int} the effect of the temperature of the internal points on the internal points
K_ib_model = Stiffness_matrix_model(idx_i_model, idx_b_model);              %% K_{int,b} the effect of the temperature of the boundary on the internal points

%% Computing A and B (implicitly)
I = eye(size(M_ii_model));
M_inv_K   = M_ii_model \ K_ii_model;   
M_inv_Kib = M_ii_model \ K_ib_model;

LHS = I + sampletime * M_inv_K;                                             %% LHS * T^{k+1} = T^k + RHS_B * Boundary_condition
RHS_B = -sampletime * M_inv_Kib;

%% Sensor measurements

sensor_pos_ratio = [0.1, 0.3, 0.5, 0.7, 0.9];                               %% Relative sensor positions 




Num_sensors = length(sensor_pos_ratio);

measured_segments = round(sensor_pos_ratio * N_Element_model);              %% Respective elements measured in model

x_sensors = x_model(measured_segments) + dx_model/2;                        %% Position of the sensor

C =  zeros(Num_sensors, N_model);                                           %% Preallocate space for sensor mapping matrix

for i = 1:Num_sensors                                                       %% Loop over all sensors

    seg_idx = measured_segments(i);                                         %% Check the segment it measures
    
    node_a = seg_idx;                                                       %% Find the corresponding nodes
    node_b = seg_idx + 1;
    


    C(i, [node_a, node_b]) = 0.5;                                           %% Both nodes contribute half to the measurement (average of the 2 nodes is measured)
end
    
%% Sensor noise
R_matrix = eye(Num_sensors);
Sensor_mean = zeros(Num_sensors, 1);                                        %% Sensor mean/bias
Y_measured = zeros(Num_sensors, length(Time_Vector));                       %% innitiate sensor measurements vector
%% Comparing process noise each step
True_T_int = True_temperature_at_Model_Nodes(idx_i_model, :);               %% Take the internal points

N_int = length(idx_i_model);                                                %% Number of internal points
W = zeros(N_int, Simulation_Steps-1);                                       %% Preallocate process noise matrix (each column is w_k)
for k = 1:Simulation_Steps-1                                                %% Calculate error during one step of the model (compared to ground truth)
    
    
    T_k_true = True_T_int(:, k);                                            %% True state at time step k
    
    
    T_k1_true = True_T_int(:, k+1);                                         %% True state at next time step k+1
    
    
    T_pred = LHS \ (T_k_true + RHS_B * Boundary_temperature_model);         %% Model prediction for time step k+1 using true state at k
    
    
    W(:, k) = T_k1_true - T_pred;                                           %% Error is process noise
    
end

w_mean = mean(W, 2);                                                        %% Remove bias from process noise
Q_est = eye(N_int) * (max(abs(w_mean)))^2;
if Forcing_known_mode
    Q_est = (sampletime)^2 * eye(N_int) * 1e-4;
end

C_int = C(:, idx_i_model);                                                  %% Compute C for the internal points only
Kalman_update_interval = sampletime;                                        %% time in seconds before each Kalman update
Steps_before_Kalman = round(Kalman_update_interval/sampletime);             %% Steps before each Kalman update
Q_factor = 1;
Q_used = Q_factor * Q_est;

b_term =  RHS_B * Boundary_temperature_model;
f_i_model = forcing_model(idx_i_model);  

Bf = LHS \ (sampletime * (M_ii_model \ eye(size(M_ii_model))));
forcing_term_model = Bf * f_i_model;    
sensor_noise_all = mvnrnd(Sensor_mean, R_matrix, Simulation_Steps * N_sim);
sensor_noise_all = reshape(sensor_noise_all', Num_sensors, Simulation_Steps, N_sim);
process_noise_all = mvnrnd(zeros(N_int, 1), Q_used, Simulation_Steps * N_sim);
process_noise_all = reshape(process_noise_all', N_int, Simulation_Steps, N_sim);
Y_true = C * True_temperature_at_Model_Nodes;

error = zeros(N_model, Simulation_Steps, N_sim);

[L,U] = lu(LHS);                                                            %% one-time factorisation
initial_error_var = 0.1;                                                    %% Variance of initial temperature error [Kelvin^2]
initial_error_covariance = eye(N_int) * initial_error_var;
for run = 1:N_sim
 T_est_full = zeros(N_model, Simulation_Steps);                             %% Preallocate full estimate

%% new Measurement
Y_measured =  Y_true + sensor_noise_all(:, :, run);

    
T_est = zeros(N_int, Simulation_Steps);                                     %% Preallotacte state estimate
P = eye(N_int) * initial_error_var;                                         %% Initial covariance change if needed

T_est(:,1) = T_0(idx_i_model) + mvnrnd(zeros(N_int,1), initial_error_covariance)'; %% Initial state

for k = 1:Simulation_Steps-1
    T_pred = U \ (L \ (T_est(:,k) + b_term));                               %% State prediction model 
    if Forcing_known_mode
    T_pred = T_pred + forcing_term_model;
    end
    
    X = U \ (L \ P);
    P_pred = (U' \ (L' \ X'))' + Q_used;
    
    y_pred = C_int * T_pred;                                                %% Measurement prediction
    
    y_meas = Y_measured(:, k+1);                                            %% Actual measurement
    
    
    Residual = y_meas - y_pred;                                             %% Compute residual
    
    
    S = C_int * P_pred * C_int' + R_matrix;                                 %% Residual covariance
    
    
    K = (P_pred * C_int') / S;                                              %% Kalman gain
    
    if mod(k, Steps_before_Kalman) == 0
        T_est(:,k+1) = T_pred + K * Residual;                               %% State update
    
   
        P = (eye(N_int) - K * C_int) * P_pred;                              %% Covariance update
    else
        T_est(:,k+1) = T_pred;
        P = P_pred;

    end
end
    for k = 1:Simulation_Steps                                              %% Add the boundary conditions to the estimate
    T_est_full(:,k) = [Boundary_temperature_model(1);
                       T_est(:,k);
                       Boundary_temperature_model(2)];
    end
%% Extracting data

error(:,:,run) = True_temperature_at_Model_Nodes - T_est_full;
end
if Forcing_known_mode
    save('data\all_errors_forcing_known.mat','error');
end
if Forcing_known_mode == false
    save('data\all_errors_forcing_unknown', 'error')
end