%% Script file to generate fracture networks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;
addpath(genpath('methods'));
format short
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Some parameters for multiple realizations:
rng(10)  % uncomment if you don't want to be able to replicate result (or simply change the seed)
nr_real = 100;  % number of realizations to generate
BASE_DIR = 'ensemble_1\\';   % folder to which the results are written
plot_results = false;  % will greatly slow down the code, but shows insight into how the fracture network is generated
post_process = true;  % almost always a good idea (will calc. intersections and remove duplicated nodes), might slow down the code
recalc_connectivity = true;  % if true will force prescribed connectivity by removing/adding connections (perhaps not geologically realistic)

% Domain parameters:
x_min = 0;
x_max = 1000;
y_min = 0;
y_max = 1000;

h  = 15;  % Discretization length (i.e., size of fracture segments)
alpha  = 3.0;  % Parameter for power-law distribution related to fracture length

% Dependent parameters (num_fracs and size_frac_set) on connectivity:
connectivity_dependency = true;

% Mean and std of all input parameters
mean_num_fracs = 165;
std_num_fracs = 40;
min_num_fracs = 50;
max_num_frac = 300;

min_conn = 0.05;
max_conn = 1.0;

mean_frac_per_set = 8;
std_frac_per_set = 2;
min_frac_set = 2;
max_frac_set = 16;

mean_small_frac_len = 40;
std_small_frac_len = 7.5;
min_small = 20;
max_small = 60;

mean_large_frac_len = 55;
std_large_frac_len = 10;
min_large = 30;
max_large = 85;

mean_angle_1 = 90;
std_angle_1 = 20;

mean_angle_2 = 0;
std_angle_2 = 20;

deviate_angles_newnetwork = true;
rand_angle_int = 10; % E.g.: 15 means random angle change between -15 and +15 degrees.
angle_deviation = [0;0];

% Perform restart algorithm for new fracture sets:
% Choose random new position (1), furthest away from 
% any existing fracture (2), or normal distribution with
% the mean of fracture furthest away (3):
random_restart_method = 1;

%% Initialization and assignment several parameters
ratio_connected = zeros(nr_real, 1);
number_of_fractures = zeros(nr_real, 1);
max_fractures_per_set = zeros(nr_real, 1);
min_length_small_fracs = zeros(nr_real, 1);
min_length_large_fracs = zeros(nr_real, 1);
angle_1_params = zeros(nr_real, 1);
angle_2_params = zeros(nr_real, 1);

Xb = x_min + h;   %x lower boundary
Xt = x_max - h;   %x upper boundary
Yb = y_min + h;   %y lower boundary
Yt = y_max - h;   %y upper boundary

if ~exist(BASE_DIR, 'dir')
   mkdir(BASE_DIR)
end

for ith_real = 1:nr_real
    ratio_connected(ith_real) = (max_conn - min_conn) * rand() + min_conn;
    
    number_of_fractures(ith_real) = round(normrnd(mean_num_fracs, std_num_fracs));
    max_fractures_per_set(ith_real) = round(normrnd(mean_frac_per_set, std_frac_per_set));
    
    if connectivity_dependency
        number_of_fractures(ith_real) = round(number_of_fractures(ith_real) * sqrt(ratio_connected(ith_real)));
        max_fractures_per_set(ith_real) = round(max_fractures_per_set(ith_real) * sqrt(ratio_connected(ith_real)));
    end
    
    if number_of_fractures(ith_real) < min_num_fracs
        number_of_fractures(ith_real) = min_num_fracs;
    elseif number_of_fractures(ith_real) > max_num_frac
        number_of_fractures(ith_real) = max_num_frac;
    end
    
    if max_fractures_per_set(ith_real) < min_frac_set
        max_fractures_per_set(ith_real) = min_frac_set;
    elseif max_fractures_per_set(ith_real) > max_frac_set
        max_fractures_per_set(ith_real) = max_frac_set;
    end
    
    min_length_small_fracs(ith_real) = round(normrnd(mean_small_frac_len, std_small_frac_len));
    if min_length_small_fracs(ith_real) < min_small
        min_length_small_fracs(ith_real) = min_small;
    elseif min_length_small_fracs(ith_real) > max_small
        min_length_small_fracs(ith_real) = max_small;
    end
    
    min_length_large_fracs(ith_real) = round(normrnd(mean_large_frac_len, std_large_frac_len));
    if min_length_large_fracs(ith_real) < min_large
        min_length_large_fracs(ith_real) = min_large;
    elseif min_length_large_fracs(ith_real) > max_large
        min_length_large_fracs(ith_real) = max_large;
    end
    
    angle_1_params(ith_real) = round(normrnd(mean_angle_1, std_angle_1));
    angle_2_params(ith_real) = round(normrnd(mean_angle_2, std_angle_2));
end

% Write parameter set to file:
fid1 = fopen([BASE_DIR 'parameter_set.txt'],'w+');
fprintf(fid1, 'Connect_Ratio \t Num_Fracs \t Max_Fracs_Per_Set \t Min_Length_Small_Fracs \t Min_Length_Large_Fracs Angle_1 \t Angle_2 \n');
for ith_param = 1:nr_real
     fprintf(fid1,[num2str(ratio_connected(ith_param)), '\t', ...
                   num2str(number_of_fractures(ith_param)), '\t', ...
                   num2str(max_fractures_per_set(ith_param)), '\t', ...
                   num2str(min_length_small_fracs(ith_param)), '\t', ...
                   num2str(min_length_large_fracs(ith_param)), '\t', ...
                   num2str(angle_1_params(ith_param)), '\t', ...
                   num2str(angle_2_params(ith_param)), '\t']);  
     fprintf(fid1,'\n'); 
end
fclose(fid1);

% Plot distribution of parameters:
num_bins = round(nr_real / 4);

figure()
subplot(3, 2, 1)
histogram(ratio_connected, num_bins)
title('connectivity')

subplot(3, 2, 2)
histogram(number_of_fractures, num_bins)
title('num fracs')

subplot(3, 2, 3)
histogram(max_fractures_per_set, num_bins)
title('size frac set')

subplot(3, 2, 4)
histogram(min_length_small_fracs, num_bins)
title('min len small')

subplot(3, 2, 5)
histogram(min_length_large_fracs, num_bins)
title('min len large')

subplot(3, 2, 6)
histogram(angle_1_params, num_bins)
hold on
histogram(angle_2_params, num_bins)
hold off
title('angle fracs')

%% Start main loop for each ensemble member (realization)
for ith_real = 1:nr_real
    % Read parameters for individual realization
    ratio_conn_networks = ratio_connected(ith_real);
    max_fracs_per_set = max_fractures_per_set(ith_real);
    f           =  number_of_fractures(ith_real); %18     %number of fractures 32
    Lb_small    =  min_length_small_fracs(ith_real);  %lower limit of the fracture length 8 --8  12
    Lt_small    =  min_length_large_fracs(ith_real); %upper limit of the fracture length 38 --30 16
    angle_1 = angle_1_params(ith_real);
    angle_2 = angle_2_params(ith_real);

    if plot_results
        figure(420);
        clf
        hold on
    end

    [act_frac_sys, current_connectivity] = ...
        generate_powerlaw_fracs(Lb_small, Lt_small, Xt, Xb, Yb, Yt, f, h, ...
            angle_deviation, angle_1, angle_2, plot_results, post_process, ...
            recalc_connectivity, max_fracs_per_set, x_max, y_max, ...
            deviate_angles_newnetwork, rand_angle_int, ratio_conn_networks, ...
            alpha, random_restart_method);

    % Write some output to a file:
    fid3 = fopen([BASE_DIR 'conn_real_' num2str(ith_real) '.txt'],'w+');
    fprintf(fid3,'%8.5f\n',current_connectivity);
    fclose(fid3);

    fid3 = fopen([BASE_DIR 'real_' num2str(ith_real) '.txt'],'w+');
    for ii = 1:size(act_frac_sys, 1)
        fprintf(fid3,'%8.5f \t %8.5f \t %8.5f \t %8.5f\n', act_frac_sys(ii, :));
    end
    fclose(fid3);
    clc; 
    if plot_results
        figure(420); 
        clf;
    end
end