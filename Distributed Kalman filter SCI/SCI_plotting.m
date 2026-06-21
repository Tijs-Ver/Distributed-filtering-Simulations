clc; clear; close all; 
rng(1);                                                                     %% Set RNG to compare different filters

%% Initilization
sampletime = 5;                                                           %% Sample time [s]
Simulation_Time = 100;                                                     %% Simulation duration [s]
Time_Vector = 0:sampletime:Simulation_Time;                                 %% Create time vector with all values
Simulation_Steps = length(Time_Vector);                                     %% Number of time steps
Forcing_known_mode = true;

Epsilon = 10^(-7);
if Forcing_known_mode
    Sample_standard_deviations = load("data\Standard_deviation_sample_known_forcing_SCI.mat").stds;    %% Loads sample standard deviations
    Mean_final_error = load("data\Mean_final_error_known_forcing_SCI.mat").Mean_final_error;
    Mean_error_over_time = load("data\Mean_error_over_time_known_forcing_SCI.mat").Mean_error_over_time;
end
if Forcing_known_mode == false
    Sample_standard_deviations = load("data\Standard_deviation_sample_unknown_forcing_SCI.mat").stds;    %% Loads sample standard deviations
    Mean_final_error = load("data\Mean_final_error_unknown_forcing_SCI.mat").Mean_final_error;
    Mean_error_over_time = load("data\Mean_error_over_time_unknown_forcing_SCI.mat").Mean_error_over_time;
end

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

%% Visualistation
%% Colorblind friendly colors
color.black  = [0 0 0]/255;
color.orange = [230 159 0]/255;
color.cyan   = [86 180 233]/255;
color.green  = [0 158 115]/255;
color.yellow = [240 228 66]/255;
color.blue   = [0 114 178]/255;
color.red    = [213 94 0]/255;
color.pink   = [204 121 167]/255;

%% Academic Plotting Style Variables
fontSizeLabel = 20;          
fontSizeLegend = 16;        
lineWidthThin = 1.5;        
lineWidthThick = 2.5;        
markerSize = 8;              
figPosition = 4*[0 0 192 144]; 
interpreterType = 'latex';   



%% Set and create folder for plots
if Forcing_known_mode
    plot_folder = 'plots_forcing_known';
else
    plot_folder = 'plots_forcing_unknown';
end

if ~exist(plot_folder,'dir')
    mkdir(plot_folder);
end

%% Heatmap: Ground truth temperature evolution
figure('Position', figPosition);
hold on; grid on; box on;
imagesc(Time_Vector, x, True_temperature);
set(gca, 'YDir', 'normal');
colormap(turbo);
clim([298 305]);
cb = colorbar;
cb.TickLabelInterpreter = 'latex';
xlabel('$\tau\,(\mathrm{s})$', 'Interpreter', 'latex');
ylabel('$x\,(\mathrm{m})$', 'Interpreter', 'latex');
ylabel(cb, '$t\,(\mathrm{K})$', 'Interpreter', 'latex');
set(gca, 'FontSize', fontSizeLabel);
set(gca, 'TickLabelInterpreter', 'latex');
hold off;
savefig(fullfile(plot_folder, 'ground_truth_heatmap.fig'));
set(gcf, 'renderer', 'Painters');
saveas(gcf, fullfile(plot_folder, 'ground_truth_heatmap.eps'), 'epsc');

%% Surface plot: ground truth temperature evolution
figure('Position', figPosition);
hold on; grid on; box on;
surf(x, Time_Vector, True_temperature');
shading interp;
colormap(turbo);
cb = colorbar;
cb.TickLabelInterpreter = 'latex';
xlabel('$x\,(\mathrm{m})$', 'Interpreter', 'latex');
ylabel('$\tau\,(\mathrm{s})$', 'Interpreter', 'latex');
zlabel('$t\,(\mathrm{K})$', 'Interpreter', 'latex');
ylabel(cb, '$t\,(\mathrm{K})$', 'Interpreter', 'latex');
set(gca, 'FontSize', fontSizeLabel);
set(gca, 'TickLabelInterpreter', 'latex');
view(-45, 30);
hold off;
savefig(fullfile(plot_folder, 'ground_truth_surface.fig'));
set(gcf, 'renderer', 'Painters');
saveas(gcf, fullfile(plot_folder, 'ground_truth_surface.eps'), 'epsc');


%% Temperature snapshots at selected times

times_to_plot = [10, 40, 70, 100];  % [s]

idx = zeros(size(times_to_plot));
for i = 1:length(times_to_plot)
    [~, idx(i)] = min(abs(Time_Vector - times_to_plot(i)));
end

%% Consistent color palette (same as thesis figures)
lineColors = {color.blue, color.green, color.orange, color.red};

figure('Position', figPosition);
hold on; grid on; box on;

for i = 1:length(idx)
    plot(x, True_temperature(:, idx(i)), ...
        'LineWidth', lineWidthThick, ...
        'Color', lineColors{i}, ...
        'DisplayName', ['$\tau = ', num2str(times_to_plot(i)), '\,\mathrm{s}$']);
end

xlabel('$x\,(\mathrm{m})$', 'Interpreter', 'latex');
ylabel('$t\,(\mathrm{K})$', 'Interpreter', 'latex');

legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', fontSizeLegend);

set(gca, 'FontSize', fontSizeLabel);
set(gca, 'TickLabelInterpreter', 'latex');

hold off;

%% Save figure (thesis format)
savefig(fullfile(plot_folder, 'ground_truth_snapshots.fig'));
set(gcf, 'renderer', 'Painters');
saveas(gcf, fullfile(plot_folder, 'ground_truth_snapshots.eps'), 'epsc');
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

sensor_var = [0.01, 0.005, 0.001, 0.005, 0.01];                                   %% Variances of each of the sensors


R_matrix = eye(Num_sensors);
Sensor_mean = zeros(Num_sensors, 1);

measurement_noise = mvnrnd(Sensor_mean, R_matrix, Simulation_Steps)';


Y_measured = zeros(Num_sensors, Simulation_Steps);

%% Final Measurement
for i = 1:Simulation_Steps
    Y_measured(:, i) = C * True_temperature_at_Model_Nodes(:, i) + measurement_noise(:, i);
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
%% Computing A and B (implicitly)
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
for k = 2:Simulation_Steps                                                %% Calculate error during one step of the model (compared to ground truth)
    
    
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
X_bound_local_Correlated = cell(N_agents, 1);
X_bound_local_Uncorrelated = cell(N_agents, 1);
X_history = cell(N_agents, 1);                                              %% Storage for estimates (all x_i_hat(k|k))
X_bound_history = cell(N_agents, 1);                                        %% Storage for all bounds (X_ii(k|k)
X_bound_history_Correlated = cell(N_agents, 1);
X_bound_history_Uncorrelated = cell(N_agents, 1);
Z_X_i = cell(N_agents, 1);                                                  %% Temporary holder of information([x_hat(k|k-1);y(k)]
x_hat_pred = cell(N_agents, 1);                                             %% Temporary holder for x_hat(k|k-1)
K_history = cell(N_agents, 1);
R_X_I_SCI = cell(N_agents, 1);
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
initial_error_covariance = eye(N_int) * initial_error_var;
Initial_temp_model = T_0_internal + mvnrnd(zeros(N_int, 1), initial_error_covariance)';
for i = 1:N_agents
    n_i = length(agent_idx{i});                                             %% Number of nodes for this agent
   
    x_hat_local{i} = Initial_temp_model(agent_idx{i}); %% initial condition
    
    X_bound_local{i} = initial_error_var * eye(n_i);                        %% initial covariance upper bound (taken to be very large just saying that we don't really have an idea about an upper bound at the start it will automatically become a better bound over time)
    X_bound_local_Uncorrelated{i} = initial_error_var * eye(n_i);
    X_bound_local_Correlated{i} = Epsilon * eye(n_i);

    X_history{i} = zeros(n_i, Simulation_Steps);                            %% Allocating space for all estimates
    X_history{i}(:, 1) = x_hat_local{i};                                    %% Storing first estimate
    X_bound_history{i} = zeros(n_i, n_i, Simulation_Steps);                 %% Allocating space for all bounds
    X_bound_history{i}(:,:,1) = X_bound_local{i};                           %% Storing initial (guess) bound
    X_bound_history_Correlated{i} = zeros(n_i, n_i, Simulation_Steps);                 %% Allocating space for all bounds
    X_bound_history_Correlated{i}(:,:,1) = X_bound_local_Correlated{i};                           %% Storing initial (guess) bound
    X_bound_history_Uncorrelated{i} = zeros(n_i, n_i, Simulation_Steps);                 %% Allocating space for all bounds
    X_bound_history_Uncorrelated{i}(:,:,1) = X_bound_local_Uncorrelated{i};                           %% Storing initial (guess) bound
    K_history{i} = zeros(n_i, n_i + 1, Simulation_Steps);
end




%% Main Distributed Estimation Loop
for k = 2:Simulation_Steps

    
    for i = 1:N_agents                                                      %% Prediction step each agent finds (x_hat(k|k-1)
  
        x_hat_pred{i} = agent_A{i,i} * x_hat_local{i} + agent_B{i} * Boundary_temperature_model; %% predict using own previous estimate
        
        
        if i > 1
            
            x_hat_pred{i} = x_hat_pred{i} + agent_A{i, i-1} * x_hat_local{i-1}; %% predict using previous estimate of left neighbour
        end
        if i < N_agents
            
            x_hat_pred{i} = x_hat_pred{i} + agent_A{i, i+1} * x_hat_local{i+1}; %% predict using previous estimate of right neighbour
        end
        if Forcing_known_mode
            x_hat_pred{i} = x_hat_pred{i} + agent_forcing_term{i};
        end
        n_i = length(agent_idx{i});                                         %% number of nodes in agent
        Z_X_i{i} = [x_hat_pred{i};                                          %% Find information vector
            Y_measured(i, k)];
        
    end
    for i = 1:N_agents                                                      %% This part is very much hard coded, it would be better if somewhere a variable that shows the communication between agents is used
        if i == 1
            Yb = cell(2, 1);
            Yb{1} = omega_b{i}{1}' * (X_bound_local_Correlated{i}\omega_b{i}{1});
            Yb{2} = omega_b{i}{2}' * (X_bound_local_Correlated{i+1}\omega_b{i}{2});
            R_X_I_SCI{i} = R_X_i{i} + C_X_i{i} * blkdiag(X_bound_local_Uncorrelated{i}, X_bound_local_Uncorrelated{i+1}) * C_X_i{i}';
            
        elseif i == N_agents
            Yb = cell(2,1);
            Yb{1} = omega_b{i}{1}' * (X_bound_local_Correlated{i-1}\omega_b{i}{1});
            Yb{2} = omega_b{i}{2}' * (X_bound_local_Correlated{i}\omega_b{i}{2});
            R_X_I_SCI{i} = R_X_i{i} + C_X_i{i} * blkdiag(X_bound_local_Uncorrelated{i-1}, X_bound_local_Uncorrelated{i}) * C_X_i{i}';
        else
            Yb = cell(3, 1);
            Yb{1} = omega_b{i}{1}' * (X_bound_local_Correlated{i-1}\omega_b{i}{1});
            Yb{2} = omega_b{i}{2}' * (X_bound_local_Correlated{i}\omega_b{i}{2});
            Yb{3} = omega_b{i}{3}' * (X_bound_local_Correlated{i+1}\omega_b{i}{3});
            R_X_I_SCI{i} = R_X_i{i} + C_X_i{i} * blkdiag(X_bound_local_Uncorrelated{i-1}, X_bound_local_Uncorrelated{i}, X_bound_local_Uncorrelated{i+1}) * C_X_i{i}';
        end
        [K, bound, ~, ~] = overlappingCI(H_X_i{i}, R_X_I_SCI{i}, C_X_i{i}, Yb, 'trace', struct('warning_on_numerical_problems', 1, 'normalization', 1)); %% Solve optimization problem , 'normalization', 1
        x_new{i} = K * Z_X_i{i};                                      %% Find estimate
        X_history{i}(:,k) = x_new{i};                                 %% Save estimate

        X_new{i} = bound;                                           %% New local bound
        X_bound_history{i}(:,:,k) = X_new{i};                       %% Save local bound

        K_history{i}(:, :, k) = K;

        X_new_Uncorrelated{i} = K * R_X_i{i} * K';
        X_bound_history_Uncorrelated{i}(:,:,k) = X_new_Uncorrelated{i};

        X_new_Correlated{i} = X_new{i} - X_new_Uncorrelated{i};
        X_bound_history_Correlated{i}(:,:,k) = X_new_Correlated{i};
        


    end
   x_hat_local = x_new;
   X_bound_local = X_new;
   X_bound_local_Uncorrelated = X_new_Uncorrelated;
   X_bound_local_Correlated = X_new_Correlated;
end

%% Post Processing and visualisation


Full_Estimate_temporary = NaN(N_int, Simulation_Steps);                       %% Initilising full estimate
Sigma3_Bound_temporary  = NaN(N_int, Simulation_Steps);                       %% Initilising full sigma 3 bound

Full_Estimate_History = NaN(N_model, Simulation_Steps);                       %% Initilising full estimate
Sigma3_Bound_History  = NaN(N_model, Simulation_Steps);                       %% Initilising full sigma 3 bound


for i = 1:N_agents
    
    
    Full_Estimate_temporary(agent_idx{i}, :) = X_history{i};                  %% obtain full estimate
    
    
    for k = 1:Simulation_Steps

        local_variance = diag(X_bound_history{i}(:, :, k));                 %% Extract local variances
        
        Sigma3_Bound_temporary(agent_idx{i}, k) = 3 * sqrt(local_variance);   %% variancse along entire rod
    end
end
for i = 1:Simulation_Steps
    Full_Estimate_History(:, i) = [Boundary_temperature_model(1);Full_Estimate_temporary(:,i);Boundary_temperature_model(2)];
    Sigma3_Bound_History(:, i) = [0;Sigma3_Bound_temporary(:,i);0];
end


x_internal_model = x_model(idx_i_model);                                    %% X positions of internal points for plotting

%% 3D plot: estimated temperature evolution (CI)

figure('Position', figPosition);
hold on; grid on; box on;

surf(x_model, Time_Vector, Full_Estimate_History');
shading interp;

colormap(turbo);

cb = colorbar;
cb.TickLabelInterpreter = 'latex';
ylabel(cb, '$t\,(\mathrm{K})$', 'Interpreter', 'latex');

xlabel('$x\,(\mathrm{m})$', 'Interpreter', 'latex');
ylabel('$\tau\,(\mathrm{s})$', 'Interpreter', 'latex');
zlabel('$t\,(\mathrm{K})$', 'Interpreter', 'latex');

set(gca, 'FontSize', fontSizeLabel);
set(gca, 'TickLabelInterpreter', 'latex');

view(-37.5, 30);

hold off;

%% Save figure (thesis format)
savefig(fullfile(plot_folder, 'ci_estimate_surface.fig'));
set(gcf, 'renderer', 'Painters');
saveas(gcf, fullfile(plot_folder, 'ci_estimate_surface.eps'), 'epsc');

Final_True     = True_temperature_at_Model_Nodes(:, end);                                        %% Final internal nodes temperature of ground truth
Final_Estimate = Full_Estimate_History(:, end);                             %% Final internal nodes temperature of estimate
Final_Error    = Final_True - Final_Estimate;                               %% Error along rod at final time
Final_3Sigma   = Sigma3_Bound_History(:, end);                              %% Bound at final time
Final_3Sigma_MC = 3 * Sample_standard_deviations(:, end);
%% Final time estimate vs ground truth

figure('Position', figPosition);
hold on; grid on; box on;

%% Ground truth (black solid)
plot(x_model, Final_True, ...
    'Color', color.black, ...
    'LineWidth', lineWidthThick, ...
    'LineStyle', '-', ...
    'DisplayName', 'Ground truth');

%% Estimate (blue dotted — consistent with Kalman-style plots)
plot(x_model, Final_Estimate, ...
    'Color', color.blue, ...
    'LineWidth', lineWidthThick, ...
    'LineStyle', ':', ...
    'DisplayName', 'Estimate');

%% Sensor positions (green dashed vertical lines)
for i = 1:length(x_sensors)
    h = xline(x_sensors(i), '--', ...
        'Color', color.green, ...
        'LineWidth', lineWidthThin);

    if i == 1
        h.DisplayName = 'Sensor positions';
    else
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end

xlabel('$x\,(\mathrm{m})$', 'Interpreter', 'latex');
ylabel('$t\,(\mathrm{K})$', 'Interpreter', 'latex');

legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', fontSizeLegend);

set(gca, 'FontSize', fontSizeLabel);
set(gca, 'TickLabelInterpreter', 'latex');

hold off;

%% Save figure (thesis format)
savefig(fullfile(plot_folder, 'final_estimate_vs_ground_truth.fig'));
set(gcf, 'renderer', 'Painters');
saveas(gcf, fullfile(plot_folder,'final_estimate_vs_ground_truth.eps'), 'epsc');
%% Final time error + bounds

figure('Position', figPosition);
hold on; grid on; box on;

%% Analytical ±3σ bounds (red dotted)
plot(x_model, Final_3Sigma + Mean_final_error, ...
    'Color', color.red, ...
    'LineWidth', lineWidthThick, ...
    'LineStyle', ':', ...
    'DisplayName', '$\pm 3\sigma$ bound (analytical)');

plot(x_model, -Final_3Sigma + Mean_final_error, ...
    'Color', color.red, ...
    'LineWidth', lineWidthThick, ...
    'LineStyle', ':', ...
    'HandleVisibility', 'off');

%% Monte Carlo bounds (blue dotted — placeholder)
plot(x_model, Final_3Sigma_MC + Mean_final_error, ...
     'Color', color.blue, ...
     'LineWidth', lineWidthThin, ...
     'LineStyle', ':', ...
     'DisplayName', '$\pm 3\sigma$ bound (Monte Carlo)');

plot(x_model, -Final_3Sigma_MC + Mean_final_error, ...
     'Color', color.blue, ...
     'LineWidth', lineWidthThin, ...
     'LineStyle', ':', ...
     'HandleVisibility', 'off');

%% Realized error (black solid)
plot(x_model, Final_Error, ...
    'Color', color.black, ...
    'LineWidth', lineWidthThick, ...
    'LineStyle', '-', ...
    'DisplayName', 'Prediction error');

xlabel('$x\,(\mathrm{m})$', 'Interpreter', 'latex');
ylabel('$t_{\mathrm{error}}\,(\mathrm{K})$', 'Interpreter', 'latex');

legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', fontSizeLegend);

set(gca, 'FontSize', fontSizeLabel);
set(gca, 'TickLabelInterpreter', 'latex');
hold off;
savefig(fullfile(plot_folder,'prediction_error_bounds.fig'));
set(gcf, 'renderer', 'Painters');
saveas(gcf, fullfile(plot_folder,'prediction_error_bounds.eps'), 'epsc');

%% Middle node error + bounds

figure('Position', figPosition);
hold on; grid on; box on;

middle_True     = True_temperature_at_Model_Nodes(round(0.5*end), :);
middle_Estimate = Full_Estimate_History(round(0.5*end), :);
middle_Error    = middle_True - middle_Estimate;
middle_3Sigma   = Sigma3_Bound_History(round(0.5*end), :);
middle_3Sigma_MC = 3 * Sample_standard_deviations(round(end/2), :);
%% Analytical ±3σ bounds (red dotted)
plot(Time_Vector, middle_3Sigma + Mean_error_over_time, ...
    'Color', color.red, ...
    'LineWidth', lineWidthThick, ...
    'LineStyle', ':', ...
    'DisplayName', '$\pm 3\sigma$ bound (analytical)');

plot(Time_Vector, -middle_3Sigma + Mean_error_over_time, ...
    'Color', color.red, ...
    'LineWidth', lineWidthThick, ...
    'LineStyle', ':', ...
    'HandleVisibility', 'off');
%% Monte Carlo bounds (blue dotted — placeholder)
plot(Time_Vector, middle_3Sigma_MC + Mean_error_over_time, ...
    'Color', color.blue, ...
    'LineWidth', lineWidthThin, ...
    'LineStyle', ':', ...
    'DisplayName', '$\pm 3\sigma$ bound (Monte Carlo)');
plot(Time_Vector, -middle_3Sigma_MC + Mean_error_over_time, ...
    'Color', color.blue, ...
    'LineWidth', lineWidthThin, ...
    'LineStyle', ':', ...
    'HandleVisibility', 'off');

%% Realized error (black solid)
plot(Time_Vector, middle_Error, ...
    'Color', color.black, ...
    'LineWidth', lineWidthThick, ...
    'LineStyle', '-', ...
    'DisplayName', 'Prediction error');

xlabel('$\tau\,(\mathrm{s})$', 'Interpreter', 'latex');
ylabel('$t_{\mathrm{error}}\,(\mathrm{K})$', 'Interpreter', 'latex');

legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', fontSizeLegend);

set(gca, 'FontSize', fontSizeLabel);
set(gca, 'TickLabelInterpreter', 'latex');
hold off;
savefig(fullfile(plot_folder,'prediction_error_time.fig'));
set(gcf, 'renderer', 'Painters');
saveas(gcf, fullfile(plot_folder,'prediction_error_time.eps'), 'epsc');
