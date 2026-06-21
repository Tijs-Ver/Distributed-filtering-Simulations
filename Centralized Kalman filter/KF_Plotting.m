clc; clear; close all; 
rng(1);                                                                     %% Set Rng to compare different filters
%% Initilization
sampletime = 5;                                                             %% Sample time [s]
Simulation_Time = 100;                                                      %% Simulation duration [s]
Time_Vector = 0:sampletime:Simulation_Time;                                 %% Create time vector with all values
Simulation_Steps = length(Time_Vector);                                     %% Number of time steps


Forcing_known_mode = false;                                                  %% True = Forcing is known to the model

if Forcing_known_mode
    Sample_standard_deviations = load("data\Standard_deviation_sample_known_forcing.mat").stds;    %% Loads sample standard deviations
    Mean_final_error = load("data\Mean_final_error_known_forcing.mat").Mean_final_error;
    Mean_error_over_time = load("data\Mean_error_over_time_known_forcing.mat").Mean_error_over_time;
end
if Forcing_known_mode == false
    Sample_standard_deviations = load("data\Standard_deviation_sample_unknown_forcing.mat").stds;    %% Loads sample standard deviations
    Mean_final_error = load("data\Mean_final_error_unknown_forcing.mat").Mean_final_error;
    Mean_error_over_time = load("data\Mean_error_over_time_unknown_forcing.mat").Mean_error_over_time;
end
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
    
    Q = power / (rod_density * Volume_element_model * specific_heat);       %% Heating rate of the entire heated element [K/s]
    
    
    Element_forcing_vector_fine = (Q * (dx / 2) * [1; 
                                             1]);                           %% [m*K/s]
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

RHS = f_i - K_ib * Boundary_temperature;                                    %% Right hand side of eq1.15 

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
%% Heatmap
figure('Position',4*[0 0 192 144]); 
hold on; grid on; box on;
set(gca,'FontSize',20);
set(gca,'TickLabelInterpreter','latex') 
imagesc(Time_Vector, x, True_temperature);
set(gca, 'YDir', 'normal');
colormap(turbo);
clim([298 305]);                                                            %% Should be changed if different settings are used (removing leads to bad contrast in some cases)
clb = colorbar;
clb.TickLabelInterpreter = 'latex';
ylabel(clb, '$t (\mathrm{K})$', 'Interpreter', 'latex');
xlabel('$\tau (\mathrm{s})$','Interpreter','latex');
ylabel('$x (\mathrm{m})$','Interpreter','latex');
hold off;
savefig(fullfile(plot_folder,'ground_truth_heatmap.fig'));
set(gcf,'renderer','Painters');
saveas(gca,fullfile(plot_folder,'ground_truth_heatmap.eps'),'epsc');


%% Surface plot
figure('Position',4*[0 0 192 144]); 
hold on; grid on; box on;
set(gca,'FontSize',20);
set(gca,'TickLabelInterpreter','latex') 
surf(x, Time_Vector, True_temperature');
shading interp; 
colormap(turbo);
view(-45, 30); 
clb = colorbar;
clb.TickLabelInterpreter = 'latex';
ylabel('$\tau (\mathrm{s})$','Interpreter','latex');
xlabel('$x (\mathrm{m})$','Interpreter','latex');
zlabel('$t (\mathrm{K})$','Interpreter','latex');
hold off;
savefig(fullfile(plot_folder,'ground_truth_surface.fig'));
set(gcf,'renderer','Painters');
saveas(gca,fullfile(plot_folder,'ground_truth_surface.eps'),'epsc');

%% Temperature Snapshots
lineColors = {color.blue, color.green, color.orange, color.red};
times_to_plot = [10, 40, 70, 100];                                          %% Chosen for 100 s simulation
idx = zeros(size(times_to_plot));
for i = 1:length(times_to_plot)
    [~, idx(i)] = min(abs(Time_Vector - times_to_plot(i)));
end

figure('Position',4*[0 0 192 144]);
hold on; grid on; box on;
set(gca,'FontSize',20);
set(gca,'TickLabelInterpreter','latex')
for i = 1:length(idx)
    plot(x, True_temperature(:, idx(i)), 'LineWidth', 2.5, ...
        'Color', lineColors{i}, ...
        'DisplayName', ['$\tau = ', num2str(times_to_plot(i)), '\,\mathrm{s}$']);
end
legend('Location','best','Interpreter','latex','FontSize', 16);
ylabel('$t (\mathrm{K})$','Interpreter','latex');
xlabel('$x (\mathrm{m})$','Interpreter','latex');
hold off;
savefig(fullfile(plot_folder,'temperature_snapshots.fig'));
set(gcf,'renderer','Painters');
saveas(gca,fullfile(plot_folder,'temperature_snapshots.eps'),'epsc');

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

R_matrix = eye(Num_sensors);                                                %% Customize for different sensor variances

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
for k = 2:Simulation_Steps                                                  %% Calculate error during one step of the model (compared to ground truth)
    
    
    T_k_true = True_T_int(:, k-1);                                          %% True state at time step k
    
    
    T_k1_true = True_T_int(:, k);                                           %% True state at next time step k+1
    
    
    T_pred(:, k) = A* T_k_true + B * Boundary_temperature_model;            %% Model prediction for time step k+1 using true state at k
    if Forcing_known_mode
        T_pred(:, k) = T_pred(:, k) + forcing_term_model;
    end
    
    W(:, k-1) = T_k1_true - T_pred(:, k);                                   %% Error is process noise
    
end

w_mean = mean(W, 2);                                                        %% Remove bias from process noise
Q_est = eye(N_int) * (max(abs(w_mean)))^2;
if Forcing_known_mode
    Q_est = (sampletime)^2 * eye(N_int) * 1e-4;
end
%% Visualisation

%% Process Noise Mean
figure('Position', 4*[0 0 192 144]); 
hold on; grid on; box on;
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex');
plot(x_model(idx_i_model), w_mean, 'LineWidth', 2.5, 'Color', color.blue);
ylabel('$w_{\mathrm{mean}} (\mathrm{K})$', 'Interpreter', 'latex');
xlabel('$x (\mathrm{m})$', 'Interpreter', 'latex');
legend({'Process noise mean'}, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 16);
hold off;
savefig(fullfile(plot_folder,'process_noise_mean.fig'));
set(gcf, 'renderer', 'Painters');
saveas(gcf, fullfile(plot_folder,'process_noise_mean.eps'), 'epsc');


Kalman_update_interval = sampletime;                                        %% time in seconds before each Kalman update
Steps_before_Kalman = round(Kalman_update_interval/sampletime);             %% Steps before each Kalman update
Q_factors = [0.01, 0.1, 1, 10, 100];
T_est_full = zeros(N_model, Simulation_Steps, length(Q_factors));           %% Preallocate full estimate

P = zeros(N_int, N_int, Simulation_Steps, length(Q_factors));
initial_error_var = 0.1;                                                    %% Variance of initial temperature error [Kelvin^2]
initial_error_covariance = eye(N_int) * initial_error_var;
initial_condition = T_0(idx_i_model) + mvnrnd(zeros(N_int,1), initial_error_covariance)';
for q = 1:length(Q_factors)
    Q_used = Q_factors(q) * Q_est;
    T_est = zeros(N_int, Simulation_Steps);                                 %% Preallotacte state estimate

    P(:,:,1,q) = eye(N_int) * initial_error_var;                            %% Initial covariance change if needed

    T_est(:,1) = initial_condition;                                         %% Initial state

    for k = 1:Simulation_Steps-1
        P_k = P(:,:,k,q);

        T_pred = A * T_est(:,k) + B * Boundary_temperature_model;           %% State prediction model
        if Forcing_known_mode
            T_pred = T_pred + forcing_term_model;
        end
    
        P_pred = A * P_k * A' + Q_used;                                     %% Covariance prediction
    
    
    
        
    
        if mod(k, Steps_before_Kalman) == 0
            y_pred = C_int * T_pred;                                        %% Measurement prediction
    
            y_meas = Y_measured(:, k+1);                                    %% Actual measurement
    
    
            Residual = y_meas - y_pred;                                     %% Compute residual
    
    
            S = C_int * P_pred * C_int' + R_matrix;                         %% Residual covariance
    
    
            K = P_pred * C_int' / S;                                        %% Kalman gain
            T_est(:,k+1) = T_pred + K * Residual;                           %% State update
    
   
            P(:,:,k+1,q) = (eye(N_int) - K * C_int) * P_pred;               %% Covariance update
        else
            T_est(:,k+1) = T_pred;
            P(:,:,k+1, q) = P_pred;

        end
    end
    for k = 1:Simulation_Steps                                              %% Add the boundary conditions to the estimate
    T_est_full(:,k, q) = [Boundary_temperature_model(1);
                       T_est(:,k);
                       Boundary_temperature_model(2)];
    end
end




%%  Visualization


%% Kalman Estimate Comparison
figure('Position', 4*[0 0 192 144]); 
hold on; grid on; box on;
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex');
plot(x_model, True_temperature_at_Model_Nodes(:,end), ...
    'Color', color.black, 'LineWidth', 2.5, 'DisplayName', 'Ground Truth');
plot(x_model, T_est_full(:,end, 3), ':', ...
    'Color', color.blue, 'LineWidth', 2.5, ...
    'DisplayName', '$t_{\mathrm{est}}$ with $\mathbf{Q} = \mathbf{Q}_{\mathrm{est}}$');
for i = 1:length(x_sensors)
    h_xline = xline(x_sensors(i), '--', 'Color', color.green, 'LineWidth', 1.5);
    if i == 1
        h_xline.DisplayName = 'Sensor positions'; 
    else
        h_xline.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$t (\mathrm{K})$', 'Interpreter', 'latex');
xlabel('$x (\mathrm{m})$', 'Interpreter', 'latex');
hold off;
savefig(fullfile(plot_folder,'kalman_comparison.fig'));
set(gcf, 'renderer', 'Painters');
saveas(gcf, fullfile(plot_folder,'kalman_comparison.eps'), 'epsc');

%% Kalman Heatmap Estimate
figure('Position', 4*[0 0 192 144]); 
hold on; grid on; box on;
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex');
imagesc(Time_Vector, x_model, T_est_full(:,:,3));
set(gca, 'YDir', 'normal');
colormap(turbo); 
clim([298 305]);                                                            %% Should be changed if different settings are used (removing leads to bad contrast in some cases)
clb = colorbar; 
clb.TickLabelInterpreter = 'latex';
ylabel(clb, '$t_{\mathrm{est}} (\mathrm{K})$', 'Interpreter', 'latex');
xlabel('$\tau (\mathrm{s})$', 'Interpreter', 'latex');
ylabel('$x (\mathrm{m})$', 'Interpreter', 'latex');
hold off;
savefig(fullfile(plot_folder,'kalman_heatmap.fig'));
set(gcf, 'renderer', 'Painters');
saveas(gcf, fullfile(plot_folder,'kalman_heatmap.eps'), 'epsc');

%% Kalman Surface Plot Estimate
figure('Position', 4*[0 0 192 144]); 
hold on; grid on; box on;
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex');
surf(x_model, Time_Vector, T_est_full(:, :, 3)');
shading interp; 
colormap(turbo); 
view(-45, 30); 
clb = colorbar; 
clb.TickLabelInterpreter = 'latex';
ylabel(clb, '$t_{\mathrm{est}} (\mathrm{K})$', 'Interpreter', 'latex');
ylabel('$\tau (\mathrm{s})$', 'Interpreter', 'latex');
xlabel('$x (\mathrm{m})$', 'Interpreter', 'latex');
zlabel('$t_{\mathrm{est}} (\mathrm{K})$', 'Interpreter', 'latex');
hold off;
savefig(fullfile(plot_folder,'kalman_surface.fig'));
set(gcf, 'renderer', 'Painters');
saveas(gcf, fullfile(plot_folder,'kalman_surface.eps'), 'epsc');

%% Compute errors
Kalman_errors = zeros(N_model, Simulation_Steps, length(Q_factors));

for i = 1: length(Q_factors)
    Kalman_errors(:,:,i) = True_temperature_at_Model_Nodes - T_est_full(:,:,i);
end

%% Extract bounds analytically
P_final = P(:, :, end, 3);                                                  %% For the Q = Q_est case
sigma_final = [0;sqrt(diag(P_final));0];                                    %% Standard deviation at each point of the rod at the end time

sigma_sample_final = Sample_standard_deviations(:,end);

upper_bound_final = 3 * sigma_final + Mean_final_error;                     %% Upper bound  
lower_bound_final = -3 * sigma_final + Mean_final_error;                    %% Lower bound

upper_bound_sample_final = 3 * sigma_sample_final + Mean_final_error;       %% Upper bound  
lower_bound_sample_final = -3 * sigma_sample_final + Mean_final_error;      %% Lower bound

P_middle_pos = squeeze(P(round(end/2), round(end/2), :, 3));                %% For the Q = Q_est case
sigma_middle = sqrt(P_middle_pos);

sigma_sample_middle = Sample_standard_deviations(round(end/2),:);

upper_bound_middle = 3 * sigma_middle + Mean_error_over_time';              %% Upper bound  
lower_bound_middle = -3 * sigma_middle + Mean_error_over_time';             %% Lower bound

upper_bound_sample_middle = 3 * sigma_sample_middle + Mean_error_over_time; %% Upper bound  
lower_bound_sample_middle = -3 * sigma_sample_middle+ Mean_error_over_time; %% Lower bound


%% Prediction Error and 3-Sigma Bounds
figure('Position', 4*[0 0 192 144]); 
hold on; grid on; box on;
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex');
plot(x_model, Kalman_errors(:, end, 3), ...
    'Color', color.black, 'LineWidth', 2.5, ...
    'DisplayName', 'Prediction error');
plot(x_model, upper_bound_final, '--', ...
    'Color', color.red, 'LineWidth', 1.5, ...
    'DisplayName', 'Analytical $\pm 3\sigma$ bounds');
plot(x_model, lower_bound_final, '--', ...
    'Color', color.red, 'LineWidth', 1.5, ...
    'HandleVisibility', 'off');
plot(x_model, upper_bound_sample_final, '--', ...
    'Color', color.blue, 'LineWidth', 1.5, ...
    'DisplayName', 'Monte Carlo $\pm 3\sigma$ bounds');
plot(x_model, lower_bound_sample_final, '--', ...
    'Color', color.blue, 'LineWidth', 1.5, ...
    'HandleVisibility', 'off');
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$t_{\mathrm{error}} (\mathrm{K})$', 'Interpreter', 'latex');
xlabel('$x (\mathrm{m})$', 'Interpreter', 'latex');
hold off;
savefig(fullfile(plot_folder,'prediction_error_bounds.fig'));
set(gcf, 'renderer', 'Painters');
saveas(gcf, fullfile(plot_folder,'prediction_error_bounds.eps'), 'epsc');

%% Temporal Prediction Error (Middle of the rod)
figure('Position', 4*[0 0 192 144]); 
hold on; grid on; box on;
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex');
plot(Time_Vector, Kalman_errors(round(end/2), :, 3), ...
    'Color', color.black, 'LineWidth', 2.5, ...
    'DisplayName', 'Prediction error');
plot(Time_Vector, upper_bound_middle, '--', ...
    'Color', color.red, 'LineWidth', 1.5, ...
    'DisplayName', 'Analytical $\pm 3\sigma$ bounds');
plot(Time_Vector, lower_bound_middle, '--', ...
    'Color', color.red, 'LineWidth', 1.5, ...
    'HandleVisibility', 'off');
plot(Time_Vector, upper_bound_sample_middle, '--', ...
    'Color', color.blue, 'LineWidth', 1.5, ...
    'DisplayName', 'Monte Carlo $\pm 3\sigma$ bounds');
plot(Time_Vector, lower_bound_sample_middle, '--', ...
    'Color', color.blue, 'LineWidth', 1.5, ...
    'HandleVisibility', 'off');
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$t_{\mathrm{error}} (\mathrm{K})$', 'Interpreter', 'latex');
xlabel('$\tau (\mathrm{s})$', 'Interpreter', 'latex');
hold off;
savefig(fullfile(plot_folder,'prediction_error_time.fig'));
set(gcf, 'renderer', 'Painters');
saveas(gcf, fullfile(plot_folder,'prediction_error_time.eps'), 'epsc');




%% Kalman Filter Error Surface
figure('Position', 4*[0 0 192 144]); 
hold on; grid on; box on;
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex');
surf(x_model, Time_Vector, Kalman_errors(:, :, 3)');
shading interp; 
colormap(turbo); 
view(-45, 30); 
clb = colorbar; 
clb.TickLabelInterpreter = 'latex';
ylabel(clb, '$t_{\mathrm{error}} (\mathrm{K})$', 'Interpreter', 'latex');
ylabel('$\tau (\mathrm{s})$', 'Interpreter', 'latex');
xlabel('$x (\mathrm{m})$', 'Interpreter', 'latex');
zlabel('$t_{\mathrm{error}} (\mathrm{K})$', 'Interpreter', 'latex');
hold off;
savefig(fullfile(plot_folder,'kalman_error_surface.fig'));
set(gcf, 'renderer', 'Painters');
saveas(gcf, fullfile(plot_folder,'kalman_error_surface.eps'), 'epsc');

%% Final Temperature Comparison (All Q values)
figure('Position', 4*[0 0 192 144]); 
hold on; grid on; box on;
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex');
q_colors = {color.orange, color.blue, color.red, color.pink, color.cyan};
plot(x_model, True_temperature_at_Model_Nodes(:,end), 'Color', color.black, ...
    'LineWidth', 2.5, 'DisplayName', 'Ground truth');
for q = 1:length(Q_factors)
    c_idx = mod(q-1, length(q_colors)) + 1;   
    plot(x_model, T_est_full(:,end,q), '--', ...
        'Color', q_colors{c_idx}, 'LineWidth', 2, ...
        'DisplayName', ['$\mathbf{Q} = ', num2str(Q_factors(q)), '\, \mathbf{Q}_{\mathrm{est}}$']);
end
for i = 1:length(x_sensors)
    h_xline = xline(x_sensors(i), '--', 'Color', color.green, 'LineWidth', 1.5);
    if i == 1
        h_xline.DisplayName = 'Sensor positions';
    else
        h_xline.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
ylabel('$t (\mathrm{K})$', 'Interpreter', 'latex');
xlabel('$x (\mathrm{m})$', 'Interpreter', 'latex');
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 14);
hold off;
savefig(fullfile(plot_folder,'comparison_Q_all_temp.fig'));
set(gcf, 'renderer', 'Painters');
saveas(gcf, fullfile(plot_folder,'comparison_Q_all_temp.eps'), 'epsc');


%% Final Error Comparison (All Q values)
figure('Position', 4*[0 0 192 144]); 
hold on; grid on; box on;
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex');
for q = 1:length(Q_factors)
    c_idx = mod(q-1, length(q_colors)) + 1;
    
    plot(x_model, Kalman_errors(:,end,q), '--', ...
        'Color', q_colors{c_idx}, 'LineWidth', 2, ...
        'DisplayName', ['$\mathbf{Q} = ', num2str(Q_factors(q)), '\, \mathbf{Q}_{\mathrm{est}}$']);
end
for i = 1:length(x_sensors)
    h_xline = xline(x_sensors(i), '--', 'Color', color.green, 'LineWidth', 1.5);
    if i == 1
        h_xline.DisplayName = 'Sensor positions';
    else
        h_xline.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
ylabel('$t_{\mathrm{error}} (\mathrm{K})$', 'Interpreter', 'latex');
xlabel('$x (\mathrm{m})$', 'Interpreter', 'latex');
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 14);
hold off;
savefig(fullfile(plot_folder,'comparison_Q_all_error.fig'));
set(gcf, 'renderer', 'Painters');
saveas(gcf, fullfile(plot_folder,'comparison_Q_all_error.eps'), 'epsc');
