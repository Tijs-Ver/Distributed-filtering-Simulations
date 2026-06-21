clc; clear; close all; 


%% Initilization
sampletime = 5;                                                           %% Sample time [s]
Simulation_Time = 100;                                                     %% Simulation duration [s]
Time_Vector = 0:sampletime:Simulation_Time;                                 %% Create time vector with all values
Simulation_Steps = length(Time_Vector);                                     %% Number of time steps
Forcing_known_mode = true;
N_sim = 100;

%% Element properties                                                       %% Example material copper
rod_density = 8960;                                                         %% Unit [kg/m^3]
Thermal_diffusivity = 1.11e-4;                                              %% Unit [m^2/s]
specific_heat = 385;                                                        %% Unit [J/(kg * K)]

%% Geometric properties
L = 1;                                                                      %% Rod length [m]
A_cross = 1e-3;                                                           %% Cross-sectional area rod (Note this should be a lot smaller then the rod length, since only 1D conduction is considered)




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
Volume_element_model = Volume_element_fine * Models_ratio;                  %% Unit [m]
%% Geometric properties model

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
    
    Q = power / (rod_density * Volume_element_model * specific_heat);        %% Heating rate of the entire heated element [K/s]
    
    
    Element_forcing_vector_fine = (Q * (dx / 2) * [1; 
                                             1]);                            %% [m*K/s]
    Element_forcing_vector_model = (Q * (dx_model / 2) * [1; 
                                             1]); 
    
    
    model_elem_idx = round(pos_ratio * N_Element_model);                    %% Corresponding element in the model (Elements 30, 50 and 70 in the example
    fine_start_elem = (model_elem_idx - 1) * round(Models_ratio) + 1;              %% First fine element part of the model element
    fine_end_elem   = model_elem_idx * round(Models_ratio);                        %% Last fine element part of the model element
    forcing_model([model_elem_idx, model_elem_idx+1]) = forcing_model([model_elem_idx, model_elem_idx+1]) + Element_forcing_vector_model;
 
    for e = fine_start_elem:fine_end_elem
        n1 = e;                                                             %% Left node of element
        n2 = e + 1;                                                         %% Right node of element
        
        %% Assembly
        forcing_fine([n1 n2]) = forcing_fine([n1 n2]) + Element_forcing_vector_fine;
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

RHS = f_i - K_ib * Boundary_temperature;                                    %% Right hand side of eq1.15 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! precalculated to save time

dT_dt = @(t, T_int) M_ii \ (RHS - K_ii * T_int);                            %% ODE 1.15 rewritten to isolate dT/dt


T0_int = T_0(idx_i);                                                        %% Initial temperature of internal nodes

[t_out, T_int_out] = ode15s(dT_dt, Time_Vector, T0_int);                    %% Call ODE15 function

True_temperature = [ones(size(t_out))*Boundary_temperature(1), T_int_out, ones(size(t_out))*Boundary_temperature(2)]'; %% Add boundary values

%% Sensor measurements
%% Interpolating ground truth to smaller grid

True_temperature_at_Model_Nodes = zeros(N_model, Simulation_Steps);         %% Preallocate sized down ground truth

for k = 1:Simulation_Steps                                                  %% Interpolate
    True_temperature_at_Model_Nodes(:, k) = interp1(x, True_temperature(:, k), x_model, 'linear')';
end

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

R_matrix = eye(Num_sensors);
Sensor_mean = zeros(Num_sensors, 1);
Y_measured = zeros(Num_sensors, Simulation_Steps, N_sim);
for n = 1:N_sim
    measurement_noise = mvnrnd(Sensor_mean, R_matrix, Simulation_Steps)';
    Y_measured(:,:, n) = C * True_temperature_at_Model_Nodes + measurement_noise;
end





%% Model
%% Initial temperature and boundary conditions model

T_0 = ones(N_model, 1) * 300;                                               %% Initial temperature rod [K]
T_0_internal = ones(N_model-2, 1) * 300;
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
f_i_model = forcing_model(idx_i_model);
%% Computing A and B 
I = eye(size(M_ii_model));
M_inv_K   = M_ii_model \ K_ii_model;   
M_inv_Kib = M_ii_model \ K_ib_model;

LHS = I + sampletime * M_inv_K;                                             %% LHS * T^{k+1} = T^k + RHS_B * Boundary_condition (+ Bf * f) 
RHS_B = -sampletime * M_inv_Kib;


%%
True_T_int = True_temperature_at_Model_Nodes(idx_i_model, :);               %% Take the internal points
N_int = length(idx_i_model);                                                %% Number of internal points
W = zeros(N_int, Simulation_Steps-1);                                       %% Preallocate process noise matrix (each column is w_k)
A = inv(LHS);                                                               %% Precomputes matrix A 
B = A * RHS_B;
Bf = A * (sampletime * (M_ii_model \ eye(size(M_ii_model))));
C_int = C(:, idx_i_model);
forcing_term_model = Bf * f_i_model;
T_pred = zeros(N_int, Simulation_Steps);
T_pred(:, 1) = T_0_internal;
for k = 2:Simulation_Steps                                                  %% Calculate error during one step of the model (compared to ground truth)
    
    
    T_k_true = True_T_int(:, k-1);                                            %% True state at time step k
    
    
    T_k1_true = True_T_int(:, k);                                         %% True state at next time step k+1
    
    
    T_pred(:, k) = A* T_k_true + B * Boundary_temperature_model;                  %% Model prediction for time step k+1 using true state at k
    if Forcing_known_mode
        T_pred(:, k) = T_pred(:, k) + forcing_term_model;
    end
    
    W(:, k-1) = T_k1_true - T_pred(:, k);                                           %% Error is process noise
    
end

w_mean = mean(W, 2);                                                        %% Remove bias from process noise
Q_est = eye(N_int) * (max(abs(w_mean)))^2;
if Forcing_known_mode
    Q_est = (sampletime)^2 * eye(N_int) * 1e-4;
end
%% Distribution

agent_idx = {1:20, 21:40, 41:60, 61:80, 81:99};                             %% Define internal nodes corresponding to each agent note that the largest node number should be equal to N_model - 2

N_agents = length(agent_idx);

                                                                            %% Pre-allocate space for local matrices
agent_A = cell(N_agents, N_agents);                                         %% Stores A_ii and A_ij
agent_B = cell(N_agents, 1);                                                %% Local boundary input effect
agent_C = cell(N_agents, 1);                                                %% Local sensor mapping
agent_Q = cell(N_agents, 1);                                                %% Local process noise covariance
agent_R = cell(N_agents, 1);                                                %% Local measurement noise covariance
agent_forcing_term = cell(N_agents, 1);


for i = 1:N_agents                                                          %% Extract matrices
    
    
    
    
    for a = 1:N_agents                                                      %% Dynamics from left side
        agent_A{i, a} = A(agent_idx{i}, agent_idx{a});
    end
    
    
    
    agent_B{i} = B(agent_idx{i}, :);                                        %% Boundary heat input ((nearly) 0 for all agents that are not on boundary)   
    
    agent_C{i} = C_int(i, agent_idx{i});                                    %% Local sensor mapping
    
    
    agent_Q{i} = Q_est(agent_idx{i}, agent_idx{i});                         %% Local covariance of process noise !!!!!!!!!!!!!!QUESTION 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    agent_R{i} = R_matrix(i, i);                                            %% local covariance of sensor noise

    agent_forcing_term{i} = forcing_term_model(agent_idx{i});
end



%% Allocating space and naming variables
x_hat_local = cell(N_agents, 1);                                            %% Temporary holder for x_i_hat(k|k)
X_bound_local = cell(N_agents, 1);                                          %% Temporary holder for X_ii(k|k)
X_history = cell(N_agents, 1);                                              %% Storage for estimates (all x_i_hat(k|k))
X_bound_history = cell(N_agents, 1);                                        %% Storage for all bounds (X_ii(k|k)
Z_X_i = cell(N_agents, 1);                                                  %% Temporary holder of information([x_hat(k|k-1);y(k)]
x_hat_pred = cell(N_agents, 1);                                             %% Temporary holder for x_hat(k|k-1)
K_history = cell(N_agents, 1);

%% Initilising constant matrices used in CI
H_X_i = cell(N_agents, 1);                                                  %% Initilisation of H_X_i
R_X_i = cell(N_agents, 1);                                                  %% Initilisation of R_X_i
C_X_i = cell(N_agents, 1);                                                  %% Initilisation of C_X_i

for i = 1:N_agents
    n_i = length(agent_idx{i});                                             %% number of nodes per agent
    H_X_i{i} = [eye(n_i);                                                   %% Computing H_X_i
            agent_C{i}];
    R_X_i{i} = blkdiag(agent_Q{i}, agent_R{i});                             %% Computing R_X_i
    if i == 1                                                               %% Computing C_X_i
        C_X_i{i} = [agent_A{i, i}, agent_A{i, i+1};
            zeros(1, n_i), zeros(1, length(agent_idx{i+1}))];
    elseif i == N_agents
        C_X_i{i} = [agent_A{i, i-1}, agent_A{i, i};
           zeros(1, length(agent_idx{i-1})), zeros(1, n_i)];
    else
        C_X_i{i} = [agent_A{i, i-1}, agent_A{i, i}, agent_A{i, i+1};
           zeros(1, length(agent_idx{i-1})), zeros(1, n_i), zeros(1, length(agent_idx{i+1}))];
    end
end

%% Preallocating and calculating pointer matrices
omega_b = cell(N_agents, 1);                                                %% Cell of cells 


for i = 1:N_agents
    if i == 1 
        n_i = length(agent_idx{i});                                             %% number of nodes per agent
        n_i_plus = length(agent_idx{i+1});                                      %% number of nodes right neighbour
        omega_bound_agent = cell(2, 1);                                     %% From omega_b * P * omega_b^T <= X_b  (Note this is for the agents on the bound so the first cell entry gives the mapping of the X_ii on the left in a 2x2 matrix
        omega_bound_agent{1} = zeros(n_i, n_i + n_i_plus);                  %% Extracts X_11 from agent 1
        omega_bound_agent{1}(:, 1:n_i) = eye(n_i);
        omega_bound_agent{2} = zeros(n_i_plus, n_i + n_i_plus);             %% Extracts X_22 from agent 1
        omega_bound_agent{2}(:, (n_i+1):end) = eye(n_i_plus);
        omega_b{i} = omega_bound_agent;
    elseif i == N_agents
        n_i_minus = length(agent_idx{i-1});                                     %% Number of nodes left neighbour
        n_i = length(agent_idx{i});                                             %% number of nodes per agent
        omega_bound_agent = cell(2, 1);                                     %% From omega_b * P * omega_b^T <= X_b  (Note this is for the agents on the bound so the first cell entry gives the mapping of the X_ii on the left in a 2x2 matrix
        omega_bound_agent{1} = zeros(n_i_minus, n_i + n_i_minus);           %% Extracts X_44 from agent 5
        omega_bound_agent{1}(:, 1:n_i_minus) = eye(n_i_minus);
        omega_bound_agent{2} = zeros(n_i, n_i + n_i_minus);                 %% Extracts X_55 from agent 5
        omega_bound_agent{2}(:, (n_i_minus+1):end) = eye(n_i);
        omega_b{i} = omega_bound_agent;
    else
        n_i_minus = length(agent_idx{i-1});                                     %% Number of nodes left neighbour
        n_i = length(agent_idx{i});                                             %% number of nodes per agent
        n_i_plus = length(agent_idx{i+1});                                      %% number of nodes right neighbour
        omega_interior_agent = cell(3, 1);                                  %% Same logic but now for extracting 3x3 matrix
        omega_interior_agent{1} = zeros(n_i_minus, n_i_minus+n_i+n_i_plus); %% Extract X_i-1,i-1 from agent i
        omega_interior_agent{1}(:, 1:n_i_minus) = eye(n_i_minus);
        omega_interior_agent{2} = zeros(n_i, n_i_minus+n_i+n_i_plus);       %% Extract X_ii from agent i
        omega_interior_agent{2}(:, (n_i_minus+1):(n_i_minus+n_i)) = eye(n_i);
        omega_interior_agent{3} = zeros(n_i_plus, n_i_minus+n_i+n_i_plus); %% Extract X_i+1,i+1 from agent i)
        omega_interior_agent{3}(:, (n_i_minus+n_i+1):end) = eye(n_i_plus);
        omega_b{i} = omega_interior_agent;
    end
end

%% Setting up initial estimate and temporary variables
initial_error_var = 0.1;                                                      %% Variance of initial temperature error [Kelvin^2]


%% Monte Carlo simulation loop
error = zeros(N_model, Simulation_Steps, N_sim);

for run = 1:N_sim

    %% Reset initial conditions (CI)
    x_hat_local = cell(N_agents, 1);
    X_bound_local = cell(N_agents, 1);
    X_history = cell(N_agents, 1);

    for i = 1:N_agents

        n_i = length(agent_idx{i});

        initial_error_covariance = eye(n_i) * initial_error_var;

        x_hat_local{i} = T_0_internal(agent_idx{i}) + ...
            mvnrnd(zeros(n_i, 1), initial_error_covariance)';

        X_bound_local{i} = initial_error_var * eye(n_i);

        X_history{i} = zeros(n_i, Simulation_Steps);
        X_history{i}(:,1) = x_hat_local{i};
    end

    %% CI LOOP
    for k = 2:Simulation_Steps

        %% Prediction
        for i = 1:N_agents

            x_hat_pred{i} = agent_A{i,i} * x_hat_local{i} + ...
                            agent_B{i} * Boundary_temperature_model;

            if i > 1
                x_hat_pred{i} = x_hat_pred{i} + agent_A{i,i-1} * x_hat_local{i-1};
            end

            if i < N_agents
                x_hat_pred{i} = x_hat_pred{i} + agent_A{i,i+1} * x_hat_local{i+1};
            end

            if Forcing_known_mode
                x_hat_pred{i} = x_hat_pred{i} + agent_forcing_term{i};
            end

            Z_X_i{i} = [x_hat_pred{i};
                        Y_measured(i,k,run)];
        end

        %% Update
        for i = 1:N_agents

            if i == 1
                Yb = cell(2,1);
                Yb{1} = omega_b{i}{1}' * (X_bound_local{i} \ omega_b{i}{1});
                Yb{2} = omega_b{i}{2}' * (X_bound_local{i+1} \ omega_b{i}{2});

            elseif i == N_agents
                Yb = cell(2,1);
                Yb{1} = omega_b{i}{1}' * (X_bound_local{i-1} \ omega_b{i}{1});
                Yb{2} = omega_b{i}{2}' * (X_bound_local{i} \ omega_b{i}{2});

            else
                Yb = cell(3,1);
                Yb{1} = omega_b{i}{1}' * (X_bound_local{i-1} \ omega_b{i}{1});
                Yb{2} = omega_b{i}{2}' * (X_bound_local{i} \ omega_b{i}{2});
                Yb{3} = omega_b{i}{3}' * (X_bound_local{i+1} \ omega_b{i}{3});
            end

            [K, bound, ~, ~] = overlappingCI( ...
                H_X_i{i}, R_X_i{i}, C_X_i{i}, Yb, ...
                'trace', struct('warning_on_numerical_problems', 1, 'normalization', 1));

            x_hat_local{i} = K * Z_X_i{i};
            X_bound_local{i} = bound;

            X_history{i}(:,k) = x_hat_local{i};
        end
    end

    %% REBUILD FULL FIELD PROPERLY (FIXED)
    T_est_full = zeros(N_model, Simulation_Steps);

    for k = 1:Simulation_Steps

        temp_state = zeros(N_int,1);

        for i = 1:N_agents
            temp_state(agent_idx{i}) = X_history{i}(:,k);
        end

        T_est_full(:,k) = [Boundary_temperature_model(1);
                           temp_state;
                           Boundary_temperature_model(2)];
    end

    %% STORE ERROR
    error(:,:,run) = True_temperature_at_Model_Nodes - T_est_full;
    fprintf('\n################ RUN %d / %d ################\n', run, N_sim);
end

%% SAVE
if Forcing_known_mode
    save('data\all_errors_forcing_known_CI.mat','error');
else
    save('data\all_errors_forcing_unknown_CI.mat','error');
end