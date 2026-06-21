clc; clear all; close all;
%% Initilization MAKE SURE THIS MATCHES THE INITILIZATION OF SCRIPT 1
Forcing_known_mode = true;

if Forcing_known_mode
    data = load('data\all_errors_forcing_known_SCI.mat');
end
if Forcing_known_mode == false
    data = load('data\all_errors_forcing_unknown_SCI.mat');
end
All_errors = data.error;               
stds = std(All_errors, 0, 3);
Mean_final_error = mean(All_errors(:, end, :), 3);
Mean_error_over_time = mean(All_errors(round(end/2), :, :), 3);
if Forcing_known_mode
    filename = 'data\Standard_deviation_sample_known_forcing_SCI';
    filename2 = 'data\Mean_final_error_known_forcing_SCI';
    filename3 = 'data\Mean_error_over_time_known_forcing_SCI';
end
if Forcing_known_mode == false
    filename = 'data\Standard_deviation_sample_unknown_forcing_SCI';
    filename2 = 'data\Mean_final_error_unknown_forcing_SCI';
    filename3 = 'data\Mean_error_over_time_unknown_forcing_SCI';
end
save(filename, 'stds');
save(filename2, 'Mean_final_error');
save(filename3, "Mean_error_over_time"); 