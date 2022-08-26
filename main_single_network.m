%% Script file to generate fracture networks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;
addpath(genpath('methods'));
format short
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Some parameters for multiple realizations:
rng(10)  % uncomment if you don't want to be able to replicate result (or simply change the seed)
BASE_DIR = 'single_network_1\\';   % folder to which the results are written
plot_results = true;  % will greatly slow down the code, but shows insight into how the fracture network is generated
post_process = true;  % almost always a good idea (will calc. intersections and remove duplicated nodes), might slow down the code
recalc_connectivity = false;  % if true will force prescribed connectivity by removing/adding connections (perhaps not geologically realistic)

% Domain parameters:
x_min = 0;
x_max = 1500;
y_min = 0;
y_max = 1500;

num_fracs = 200;
angle_1 = 90;
angle_2 = 0;
connectivity = 0.5;

fracs_per_set = 12;
min_length_small_fracs = 40;
min_length_large_fracs = 55;

deviate_angles_newnetwork = true;
rand_angle_int = 5; % E.g.: 15 means random angle change between -15 and +15 degrees.
angle_deviation = [0;0];

% Perform restart algorithm for new fracture sets:
% Choose random new position (1), furthest away from 
% any existing fracture (2), or normal distribution with
% the mean of fracture furthest away (3):
random_restart_method = 2;

h  = 15;  % Discretization length (i.e., size of fracture segments)
alpha  = 3.0;  % Parameter for power-law distribution related to fracture length

%% Initialization and assignment several parameters
Xb = x_min + h;   %x lower boundary
Xt = x_max - h;   %x upper boundary
Yb = y_min + h;   %y lower boundary
Yt = y_max - h;   %y upper boundary

if ~exist(BASE_DIR, 'dir')
   mkdir(BASE_DIR)
end

% Write parameter set to file:
fid1 = fopen([BASE_DIR 'parameter_set.txt'],'w+');
fprintf(fid1, 'Connect_Ratio \t Num_Fracs \t Max_Fracs_Per_Set \t Min_Length_Small_Fracs \t Min_Length_Large_Fracs Angle_1 \t Angle_2 \n');
fprintf(fid1,[num2str(connectivity), '\t', ...
           num2str(num_fracs), '\t', ...
           num2str(fracs_per_set), '\t', ...
           num2str(min_length_small_fracs), '\t', ...
           num2str(min_length_large_fracs), '\t', ...
           num2str(angle_1), '\t', ...
           num2str(angle_2), '\t']);  
fprintf(fid1,'\n'); 
fclose(fid1);

%% Start main loop for generating fractures
if plot_results
    figure(420);
    clf
    hold on
end

[act_frac_sys, current_connectivity] = ...
    generate_powerlaw_fracs(min_length_small_fracs, min_length_large_fracs, Xt, Xb, Yb, Yt, ...
        num_fracs, h, angle_deviation, ...
        angle_1, angle_2, plot_results, post_process, recalc_connectivity, ...
        fracs_per_set, x_max, y_max, deviate_angles_newnetwork, ...
        rand_angle_int, connectivity, alpha, random_restart_method);

% Write some output to a file:
fid3 = fopen([BASE_DIR 'connectivity.txt'],'w+');
fprintf(fid3,'%8.5f\n',current_connectivity);
fclose(fid3);

fid3 = fopen([BASE_DIR 'frac_network_nodes.txt'],'w+');
for ii = 1:size(act_frac_sys, 1)
    fprintf(fid3,'%8.5f \t %8.5f \t %8.5f \t %8.5f\n', act_frac_sys(ii, :));
end
fclose(fid3);
clc; 