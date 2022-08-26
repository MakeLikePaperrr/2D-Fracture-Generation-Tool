% MIT License
% 
% Copyright (c) 2022 Stephan de Hoop
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
function [act_frac_sys, current_connectivity] = ...
    generate_powerlaw_fracs(len_small_fracs, len_large_fracs, Xt, Xb, Yb, Yt, ...
    f, h, angleDeviation, angle_1, angle_2, plot_results, post_process, ...
    recalc_connectivity, max_fracs_per_set, x_max, y_max, ...
    deviateAnglesOfNewNetwork, randAngleInterval, ratioOfConnectedNetworks, ...
    alpha, randomRestartMethod) 
failed_create_line_restart_tolerance = 50;
failed_create_points_restart_tolerance = 50;
failed_create_line_end_tolerance = 1000;
failed_create_points_end_tolerance = 1000;
rand_int = NaN;
min_num_segm_small = ceil(len_small_fracs / h) + 1;
max_num_segm_small = floor( (Xt - Xb) / h) - 1;
min_num_segm_large = ceil(len_large_fracs / h) + 1;
max_num_segm_large = floor( (Xt - Xb) / h) - 1;

nPointOutside=0;
nFailedToCreate=0;
restartElsewhere=false;
domainIsFull=false;

threshold_conn = 0.05;  % in fraction, not %
dist_factor = 1;
nx_restart = 10;
ny_restart = 10;
x_coord_for_restart = repmat(linspace(Xb + len_large_fracs*dist_factor, ...
    Xt - len_large_fracs*dist_factor, nx_restart), [ny_restart, 1]);
y_coord_for_restart = repmat(linspace(Yb + len_large_fracs*dist_factor, ...
    Yt - len_large_fracs*dist_factor, ny_restart)', [1, nx_restart]);

num_segm_frac_1 = zeros(f, 1);
num_segm_frac_2 = zeros(f, 1);
frac_set_vec = zeros(f, 1);
frac_set = 1;

%To assign different colors to eac fracture until 64 tones
if f<inf
    colorjet = jet;
    coll_def_size = colorjet(floor(linspace(1,64,f)),:);
else
    colorjet = jet;
    coll_def_size = colorjet;
end

current_frac_per_set = 1;
n = 1;
while (n <= f) && (~domainIsFull)
    if any(n == [round(1 * f/5), round(2 * f/5), round(3 * f/5), round(4 * f/5), f])
        fprintf('Current progress: %f \n',round(n/f * 100))
    end
    
    [mult_h, fraq_length(n)] = get_fracture_length_powerlaw( len_small_fracs, alpha, h );
    if mult_h < min_num_segm_large
        mult_h = min_num_segm_large;
        fraq_length(n) = mult_h * h;
        
    elseif mult_h > max_num_segm_large
        mult_h = max_num_segm_large;
        fraq_length(n) = mult_h * h;
    end

    [mult2_h, fraq_length2(n)] = get_fracture_length_powerlaw( len_large_fracs, alpha, h);
    if mult2_h < min_num_segm_small
        mult2_h = min_num_segm_small;
        fraq_length2(n) = mult2_h * h;
        
    elseif mult2_h > max_num_segm_small
        mult2_h = max_num_segm_small;
        fraq_length2(n) = mult2_h * h;
    end
   
    if and(mult_h >= min_num_segm_large, mult2_h >= min_num_segm_small)
        angle(n)  = angle_1+angleDeviation(1);
        angle2(n) = angle_2+angleDeviation(2);
        
        if n==1
            % Starting point of the first line of first set
            xb = Xb + (Xt-Xb)*rand;
            yb = Yb + (Yt-Yb)*rand;
            xb_new = xb;
            yb_new = yb;
        else
            if restartElsewhere
                restartElsewhere=false;
            else
                xb=[];
                yb=[];
                countLoops=0;
                while (isempty(xb) && isempty(yb) && (countLoops <=100))
                    rand_int = randi([2 length(x_disc2{n-1})-1],1);
                    xb = setdiff( x_disc2{n-1}(rand_int) , inters_x{n-1} );% = x_dim; %all the points (x) for the first segmented line of each set
                    yb = setdiff( y_disc2{n-1}(rand_int) , inters_y{n-1} );% = y_dim; %all the points (y) for the first segmented line of each set
                    countLoops=countLoops+1;
                end
            end
 
            if isempty(xb)
                xb = inters_x{n-1}; 
            end
            
            if isempty(yb)
                yb = inters_y{n-1}; 
            end
            
            xb_new = Xb - x_max;
            yb_new = Yb - y_max;
            
            countLoops=0;
            while ((xb_new < Xb) || (xb_new > Xt) || (yb_new < Yb) || (yb_new > Yt))  && (countLoops <=40)
                if countLoops <= 30
                    [xb_new,yb_new] = discretize_frac_line(xb,yb,-h,mult_h,angle(n),1);
                else
                    [xb_new,yb_new] = discretize_frac_line(xb,yb,h,mult_h,angle(n),1);
                end
                
                % Crossing of next pair with previous one
                rand_int = randi([1 length(xb_new)],1);
                xb_new = xb_new(rand_int);
                yb_new = yb_new(rand_int);
                countLoops=countLoops+1;
            end
            
        end
        
        % Generation of the first segmented line
        % Failure if xb or yb is empty
        if ~and( isempty(xb) , isempty(yb) )
            
            [x_dim,y_dim] = discretize_frac_line(xb_new,yb_new,h,mult_h,angle(n),1);
            
            % Sort from x_dim(1) < x_dim(2) < ...
            if x_dim(1) > x_dim(end)
                x_dim = flip(x_dim);
                y_dim = flip(y_dim);
            end
            
            % Check if x_dim and y_dim are fully in domain: 
            min_id = length(x_dim);
            [~, id_x] = find(x_dim > Xt);
            [~, id_y] = find(y_dim > Yt);
            if ~isempty(id_x)
                min_id = id_x(1) - 1;
            end
            if ~isempty(id_y)
                if id_y(1) < min_id
                    min_id = id_y(1) - 1;
                end
            end
            
            % Only take the part of the fracture that is inside the domain:
            if (min_id - 1) >= min_num_segm_small
                x_dim = x_dim(1:min_id);
                y_dim = y_dim(1:min_id);
                mult_h = min_id - 1;
            end
            
            fraq_length(n) = mult_h * h;
            
            x_disc1{n} = x_dim; %all the points (x) for the first segmented line of each set
            y_disc1{n} = y_dim; %all the points (y) for the first segmented line of each set
            
            rand_int_old = rand_int;
            rand_int = randi([2 mult_h],1);
            
            counter = 0;
            while (rand_int == rand_int_old) && (counter < 60)
                rand_int = randi([2 mult_h],1);
                counter = counter + 1;
            end
            
            if counter == 60
                continue
            end
            
            [x_dim, y_dim] = discretize_frac_line(x_disc1{n}(rand_int),y_disc1{n}(rand_int),h,mult2_h,angle2(n),0);
            
            % Sort from x_dim(1) < x_dim(2) < ...
            if x_dim(1) > x_dim(end)
                x_dim = flip(x_dim);
                y_dim = flip(y_dim);
            end
            
            % Check if x_dim and y_dim are fully in domain: 
            min_id = length(x_dim);
            [~, id_x] = find(x_dim > Xt);
            [~, id_y] = find(y_dim > Yt);
            if ~isempty(id_x)
                min_id = id_x(1) - 1;
            end
            if ~isempty(id_y)
                if id_y(1) < min_id
                    min_id = id_y(1) - 1;
                end
            end
            
            max_id = 1;
            [~, id_x] = find(x_dim < Xb);
            [~, id_y] = find(y_dim < Yb);
            if ~isempty(id_x)
                max_id = id_x(end) + 1;
            end
            if ~isempty(id_y)
                if id_y(1) > min_id
                    max_id = id_y(1) + 1;
                end
            end
            
            % Only take the part of the fracture that is inside the domain:
            if (min_id - max_id) >= min_num_segm_large
                x_dim = x_dim(max_id:min_id);
                y_dim = y_dim(max_id:min_id);
                mult2_h = min_id - max_id;
            end
            
            fraq_length2(n) = mult2_h * h;
            
            x_disc2{n} = x_dim;
            y_disc2{n} = y_dim;
            
            % Finally this is the list of all the coordinates of each line
            xb = x_disc1{n}(1);
            yb = y_disc1{n}(1);
            xt = x_disc1{n}(end); % resulting from fucntion for "first segmented line" saved on x_disc1{n} and the end position
            yt = y_disc1{n}(end);
            xb2 = x_disc2{n}(1); % resulting form intersection between first and second line saved on X_disc2{n} first position
            yb2 = y_disc2{n}(1);
            xt2 = x_disc2{n}(end); % resulting form function for "second segmented line" saved on X_disc2{n} end position
            yt2 = y_disc2{n}(end);
            
            inters_x{n} = x_disc1{n}(rand_int);
            inters_y{n} = y_disc1{n}(rand_int);
            
            % Cond: To restrict the points into the domain
            if xb>Xb && xb<Xt && yb>Yb && yb<Yt && xt>Xb && xt<Xt && yt>Yb && yt<Yt ...
                    &&  xb2>Xb && xb2<Xt && yb2>Yb && yb2<Yt && xt2>Xb && xt2<Xt && yt2>Yb && yt2<Yt
                % Plotting the lines
                if plot_results
                    hold on
                    grid on
                end
                
                % Cond: Check that fractures dont intersect each other and that follow 'h' spacing
                % For the first fracture n=1
                check = zeros(1,n-1);
                if n==1
                    check(1) = 1;
                end
                
                % For the second fracture on
                b = 0;
                for i=1:n-1
                    inter1 = intersection(x_disc1{i}(1),y_disc1{i}(1),x_disc1{i}(end),y_disc1{i}(end), xb2, yb2, xt2, yt2);
                    inter2 = intersection(x_disc2{i}(1),y_disc2{i}(1),x_disc2{i}(end),y_disc2{i}(end), xb_new, yb_new, xt, yt);

                    if  and(isempty(inter1), isempty(inter2)) %and(isempty(inter1), isempty(inter2))
                        check(i) = 1;
                    end
                end
                
                if  ~isempty(check) && all(check==1) %&& ~isempty(check)&& all(check==1) && b==0
                    
                    if plot_results
                        hold on
                        grid on

                        plot([x_disc1{n}(1) x_disc1{n}(end)],[y_disc1{n}(1) y_disc1{n}(end)],'Color',coll_def_size(n,:));
                        plot([x_disc2{n}(1) x_disc2{n}(end)],[y_disc2{n}(1) y_disc2{n}(end)],'Color',coll_def_size(n,:));

                        %Plotting discretized lines
                        hold on
                        plot(x_disc1{n},y_disc1{n},'.','Color', coll_def_size(n,:));
                        hold on
                        plot(x_disc2{n},y_disc2{n},'.','Color', coll_def_size(n,:));
                        drawnow;
                        
                        xlim([Xb Xt]);
                        ylim([Yb Yt]);
                    end

                    nPointOutside = 0;
                    nFailedToCreate = 0;
                    current_frac_per_set = current_frac_per_set + 1;
                    num_segm_frac_1(n) = size(x_disc1{n}, 2) - 1;
                    num_segm_frac_2(n) = size(x_disc2{n}, 2) - 1;
                    frac_set_vec(n) = frac_set;
                    n = n+1;
                else
                    nFailedToCreate = nFailedToCreate + 1;
                    
                    if nFailedToCreate > failed_create_line_end_tolerance
                        domainIsFull = true;
                        f = n-1;
                        x_disc1{f+1} = {};
                        x_disc2{f+1} = {};
                        y_disc1{f+1} = {};
                        y_disc2{f+1} = {};
                        inters_x{f+1} = {};
                        inters_y{f+1} = {};
                    elseif nFailedToCreate > failed_create_line_restart_tolerance
                        restartElsewhere = true;
                        if frac_set_vec(n-1) == frac_set
                            frac_set = frac_set + 1;
                        end
                    end
                end
                
                if current_frac_per_set >= max_fracs_per_set
                    current_frac_per_set = 0;
                    restartElsewhere = true;
                    frac_set = frac_set + 1;
                end
            else
                nPointOutside = nPointOutside + 1;
                
                if nPointOutside > failed_create_points_end_tolerance
                    domainIsFull = true;
                    f = n-1;
                    x_disc1{f+1} = {};
                    x_disc2{f+1} = {};
                    y_disc1{f+1} = {};
                    y_disc2{f+1} = {};
                    inters_x{f+1} = {};
                    inters_y{f+1} = {};
                elseif nPointOutside > failed_create_points_restart_tolerance
                    restartElsewhere = true;
                    if frac_set_vec(n-1) == frac_set
                        frac_set = frac_set + 1;
                    end
                end
                
                if current_frac_per_set >= max_fracs_per_set
                    current_frac_per_set = 0;
                    restartElsewhere = true;
                    if frac_set_vec(n-1) == frac_set
                        frac_set = frac_set + 1;
                    end
                end
            end
        end
        % If restart flag is true:
        if restartElsewhere == true
           % Perform restart algorithm:
           % Choose random new position (1), furthest away from 
           % any existing fracture (2), or normal distribution with
           % the mean of fracture furthest away (3): 
           if randomRestartMethod == 1
               xb = rand()*Xt; 
               yb = rand()*Yt; 
               
           elseif randomRestartMethod == 2 || randomRestartMethod == 3
               % Calculate point furthest away from fracture:
               % Random array with 10 x 10 points:
               dist_L_1 = 0;
               dist_R_1 = 0;
               
               dist_L_2 = 0;
               dist_R_2 = 0;
               
               dist_old = 0;
               max_dist = 0;
                
               x_coord_1_1 = zeros(n-1, 1);
               y_coord_1_1 = zeros(n-1, 1);
               x_coord_2_1 = zeros(n-1, 1);
               y_coord_2_1 = zeros(n-1, 1);
               
               x_coord_1_2 = zeros(n-1, 1);
               y_coord_1_2 = zeros(n-1, 1);
               x_coord_2_2 = zeros(n-1, 1);
               y_coord_2_2 = zeros(n-1, 1);
               
               for ith_frac = 1:n-1
                   x_coord_1_1(ith_frac) = x_disc1{ith_frac}(1);
                   y_coord_1_1(ith_frac) = y_disc1{ith_frac}(end);
                   x_coord_2_1(ith_frac) = x_disc1{ith_frac}(1);
                   y_coord_2_1(ith_frac) = y_disc1{ith_frac}(end);
                   
                   x_coord_1_2(ith_frac) = x_disc2{ith_frac}(1);
                   y_coord_1_2(ith_frac) = y_disc2{ith_frac}(end);
                   x_coord_2_2(ith_frac) = x_disc2{ith_frac}(1);
                   y_coord_2_2(ith_frac) = y_disc2{ith_frac}(end);
               end
               for ii = 1:nx_restart
                   for jj = 1:ny_restart     
                        dist_L_1 = min(sqrt( (x_coord_1_1 - x_coord_for_restart(ii, jj)).^2 + (y_coord_1_1 - y_coord_for_restart(ii, jj)).^2 ));
                        dist_R_1 = min(sqrt( (x_coord_2_1 - x_coord_for_restart(ii, jj)).^2 + (y_coord_2_1 - y_coord_for_restart(ii, jj)).^2 ));
                        
                        dist_L_2 = min(sqrt( (x_coord_1_2 - x_coord_for_restart(ii, jj)).^2 + (y_coord_1_2 - y_coord_for_restart(ii, jj)).^2 ));
                        dist_R_2 = min(sqrt( (x_coord_2_2 - x_coord_for_restart(ii, jj)).^2 + (y_coord_2_2 - y_coord_for_restart(ii, jj)).^2 ));
                        
                        min_dist = min([dist_L_1, dist_R_1, dist_L_2, dist_R_2]);
                        
                        if min_dist > dist_old
                            dist_old = min_dist;
                            x_id = ii;
                            y_id = jj;
                        end
                   end          
               end
               
               if randomRestartMethod == 2
                   xb = x_coord_for_restart(x_id, y_id); 
                   yb = y_coord_for_restart(x_id, y_id);
               end
               
               if randomRestartMethod == 3
                   xb = normrnd(x_coord_for_restart(x_id, y_id), len_large_fracs); 
                   yb = normrnd(y_coord_for_restart(x_id, y_id), len_large_fracs);
               end
                
           end

           if deviateAnglesOfNewNetwork
               angleDeviation=[randAngleInterval*(2*rand()-1);randAngleInterval*(2*rand()-1)];
           end

           xb_new = xb;
           yb_new = yb;
        end
    end
end

%% Extract act_frac_sys array (of size f_1)
tot_num_segm = sum(num_segm_frac_1) + sum(num_segm_frac_2);
act_frac_sys = zeros(tot_num_segm, 4);
act_frac_set = frac_set_vec;
frac_set_vec = zeros(tot_num_segm, 1);

for ith_segm = 1:f
    % Local to Global coordinates:
    add_1 = sum(num_segm_frac_1(1:(ith_segm - 1)));
    
    % Store fracture 1:
    for ii = 1:num_segm_frac_1(ith_segm)
        act_frac_sys(add_1 + ii, :) = [x_disc1{ith_segm}(ii), ...
                                       y_disc1{ith_segm}(ii), ...
                                       x_disc1{ith_segm}(ii + 1), ...
                                       y_disc1{ith_segm}(ii + 1)];
        frac_set_vec(add_1 + ii) = act_frac_set(ith_segm);
    end
    
    % Local to Global coordinates:
    add_2 = sum(num_segm_frac_1) + sum(num_segm_frac_2(1:(ith_segm - 1)));
    
    % Store fracture 2:
    for jj = 1:num_segm_frac_2(ith_segm)
        act_frac_sys(add_2 + jj, :) = [x_disc2{ith_segm}(jj), ...
                                       y_disc2{ith_segm}(jj), ...
                                       x_disc2{ith_segm}(jj + 1), ...
                                       y_disc2{ith_segm}(jj + 1)];
        frac_set_vec(add_2 + jj) = act_frac_set(ith_segm);
    end
end

if post_process == true || recalc_connectivity == true
% Find intersection for each segment and add extra segments (by splitting):
tolerance_intersect = 1e-4;
tolerance_zero = 1e-5;
decimals = 5;
tol_angle = 7.5;

act_frac_sys = round(act_frac_sys * 10^decimals) * 10^(-decimals);
store_act_frac_sys = act_frac_sys;
store_frac_set_vec = frac_set_vec;

% First straighten fractures:
straighten_fracs_script;

act_frac_sys = round(act_frac_sys * 10^decimals) * 10^(-decimals);

% Make sure segments are sorted from x_min to x_max:
act_frac_sys_dummy = act_frac_sys;
indices = act_frac_sys_dummy(:, 1) > act_frac_sys_dummy(:, 3);
act_frac_sys(indices, 1:2) = act_frac_sys_dummy(indices, 3:4);
act_frac_sys(indices, 3:4) = act_frac_sys_dummy(indices, 1:2);

[act_frac_sys, id_frac_set, ~] = unique(act_frac_sys, 'rows');
frac_set_vec = frac_set_vec(id_frac_set);

calc_intersections;

act_frac_sys = round(act_frac_sys * 10^decimals) * 10^(-decimals);

% Make sure segments are sorted from x_min to x_max:
act_frac_sys_dummy = act_frac_sys;
indices = act_frac_sys_dummy(:, 1) > act_frac_sys_dummy(:, 3);
act_frac_sys(indices, 1:2) = act_frac_sys_dummy(indices, 3:4);
act_frac_sys(indices, 3:4) = act_frac_sys_dummy(indices, 1:2);

[act_frac_sys, id_frac_set, ~] = unique(act_frac_sys, 'rows');
frac_set_vec = frac_set_vec(id_frac_set);
end

% Find unique nodes:
unique_nodes = [act_frac_sys(:, [1, 2]); act_frac_sys(:, [3, 4])];
[unique_nodes, ~, ~] = unique(unique_nodes, 'rows');

% Store for each node to which segments it belongs and what its degree is:
num_unq_nodes = length(unique_nodes);
num_main_segm = size(act_frac_sys, 1);
degree_nodes = zeros(num_unq_nodes, 1);
incidence_matrix = zeros(num_unq_nodes, num_main_segm);

for ii = 1:num_unq_nodes
    ids_1 = find( ((act_frac_sys(:, 1) - unique_nodes(ii, 1)).^2 + ...
                   (act_frac_sys(:, 2) - unique_nodes(ii, 2)).^2).^(1/2) < 1e-3);
    ids_2 = find( ((act_frac_sys(:, 3) - unique_nodes(ii, 1)).^2 + ...
                   (act_frac_sys(:, 4) - unique_nodes(ii, 2)).^2).^(1/2) < 1e-3);
    degree_nodes(ii) = size(ids_1, 1) +  size(ids_2, 1);
    incidence_matrix(ii, union(ids_1, ids_2)) = 1;
end

% Calculate Adjacency matrix:
degree_matrix = diag(sum(incidence_matrix, 2));

adjacency_matrix = incidence_matrix * incidence_matrix' - degree_matrix;

laplacian_matrix = degree_matrix - adjacency_matrix;

component_basis = null(laplacian_matrix, 'r');

num_comp = size(component_basis, 2);

if num_comp == 0
    num_comp = 1;
    component_basis = ones(size(laplacian_matrix, 1), 1);
end

% Store each component:
components_network_nodes = zeros(num_unq_nodes, num_comp);

% Component segments:
components_network_segms = zeros(num_main_segm, num_comp);
component_segm_tot_length = zeros(num_comp, 1);
component_num_segm = zeros(num_comp, 1);

for ii = 1:num_comp
    ids = find(abs(component_basis(:, ii) - 1) < tolerance_zero);
    components_network_nodes(ids, :) = ii;
end

for ii = 1:num_comp
    comp_segm_list = [];
    ids = find(components_network_nodes(:, ii) == ii);
    comp_segm_mat = incidence_matrix(ids, :);
    for jj = 1:size(comp_segm_mat, 1)
        comp_segm_list = [comp_segm_list, find(comp_segm_mat(jj, :))];
    end
    comp_segm_list = unique(comp_segm_list);
    components_network_segms(comp_segm_list, ii) = ii;
    
    % Calculate component total length:
    component_num_segm(ii) = length(comp_segm_list);
    component_segm_tot_length(ii) = sum( ((act_frac_sys(comp_segm_list, 1) - act_frac_sys(comp_segm_list, 3)).^2 + ...
                                          (act_frac_sys(comp_segm_list, 2) - act_frac_sys(comp_segm_list, 4)).^2).^(1/2) );
end

if plot_results
    % Plot components of graph:
    colorjet = jet;
    coll_def_size = colorjet(floor(linspace(1, 64, num_comp)),:);

    figure()
    hold on

    for ii = 1:num_comp
        ids = find(components_network_segms(:, ii) == ii);
        plot(act_frac_sys(ids, [1, 3])', act_frac_sys(ids, [2, 4])', 'color', coll_def_size(ii, :))    
    end
    hold off
end

% Find closest distances to 
[max_size, max_id] = max(component_segm_tot_length);
max_connected_ratio = max_size / sum(component_segm_tot_length);
current_connectivity = max_connected_ratio;

if recalc_connectivity == true
components_dist_ids = zeros(num_comp, num_comp, 2);  % node from ii and jj
components_dist_val = zeros(num_comp, num_comp);

for ii = 1:num_comp
    for jj = (ii + 1):num_comp
        % Choose smallest cluster as embedded loop:
        min_dist = inf;
        node_id_jj = inf;
        node_id_ii = inf;
        
        ids_ii = find(components_network_nodes(:, ii) == ii)';
        ids_jj = find(components_network_nodes(:, jj) == jj)';
        
        if component_num_segm(ii) > component_num_segm(jj)
            % Loop over jj:
            for kk = ids_jj
                [new_dist, node] = min( ( (unique_nodes(ids_ii, 1) - unique_nodes(kk, 1)).^2 + ...
                                          (unique_nodes(ids_ii, 2) - unique_nodes(kk, 2)).^2 ).^(1/2)  );
                if new_dist < min_dist
                    % Store new distance value and which node:
                    min_dist = new_dist;
                    node_id_jj = kk;
                    
                    % Find global ID for ii:
                    node_id_ii = ids_ii(node);
                end
            end
        else
            % Loop over ii:
            for kk = ids_ii
                [new_dist, node] = min( ( (unique_nodes(ids_jj, 1) - unique_nodes(kk, 1)).^2 + ...
                                          (unique_nodes(ids_jj, 2) - unique_nodes(kk, 2)).^2 ).^(1/2)  );
                if new_dist < min_dist
                    % Store new distance value and which node:
                    min_dist = new_dist;
                    node_id_ii = kk;
                    
                    % Find global ID for jj:
                    node_id_jj = ids_jj(node);
                    
                end 
            end
        end
        
        components_dist_ids(ii, jj, :) = [node_id_ii, node_id_jj];
        components_dist_val(ii, jj) = min_dist;
    end
end

% Determine which components to merge:
% Start from largest component and start adding until connectivity reached:
new_segm_to_add = [];
comp_added = [max_id];

comp_dist_val_tot = eye(num_comp) * 999999 + components_dist_val + components_dist_val';
comp_dist_ids_tot = components_dist_ids;
comp_dist_ids_tot(:, :, 1) = comp_dist_ids_tot(:, :, 1) + comp_dist_ids_tot(:, :, 1)';
comp_dist_ids_tot(:, :, 2) = comp_dist_ids_tot(:, :, 2) + comp_dist_ids_tot(:, :, 2)';

while (current_connectivity + threshold_conn) < ratioOfConnectedNetworks
    % At each iteration add component to graph:
    min_dist = inf;
    id = inf;
    
    comp_check = 1:num_comp;
    comp_check(comp_added) = [];
    
    % Find closest sub-graph to maximum graph:
    for ii = comp_check
        new_dist = comp_dist_val_tot(max_id, ii);
        
        if new_dist < min_dist
            min_dist = new_dist;
            id = ii;
        end
    end
    
    % Update added component to maximum component:
    comp_added = [comp_added; id];
    
    % Add id to max_id component:
    nodes_to_add = reshape(comp_dist_ids_tot(max_id, id, :), 1, 2);
    new_segm_to_add = [new_segm_to_add; unique_nodes(nodes_to_add(1), :), ...
                                        unique_nodes(nodes_to_add(2), :)];
    
    % Update current connectivity:
    current_connectivity = current_connectivity + component_segm_tot_length(id) / sum(component_segm_tot_length);
    
    % Change distance based on newly added component:
    comp_dist_val_tot(max_id, id) = 999999;
    comp_dist_val_tot(id, max_id) = 999999;
    
    % Choose smallest distance between maximum and added component and
    % store that in maximum component distance and ids matrix:
    comp_check = 1:num_comp;
    comp_check(comp_added) = [];
    
    for ii = comp_check
        if comp_dist_val_tot(max_id, ii) > comp_dist_val_tot(id, ii)
            % Update distance:
            comp_dist_val_tot(max_id, ii) = comp_dist_val_tot(id, ii);
            comp_dist_val_tot(ii, max_id) = comp_dist_val_tot(id, ii);

            % Update node pairs:
            comp_dist_ids_tot(max_id, ii, :) = comp_dist_ids_tot(id, ii, :);
            comp_dist_ids_tot(ii, max_id, :) = comp_dist_ids_tot(id, ii, :);
        end
    end
end

% Check if connectivit is far above target connectivity:
ratio_connected_threshold = current_connectivity / ratioOfConnectedNetworks;
reduction_in_length = 1 / ratio_connected_threshold;
  
% Find length of fracture sets:
num_sets = max(frac_set_vec);
length_frac_sets = zeros(num_sets, 1);

for ii = 1:num_sets
    length_frac_sets(ii) = sum( ((act_frac_sys(frac_set_vec==ii, 1) - act_frac_sys(frac_set_vec==ii, 3)).^2 + ...
                                 (act_frac_sys(frac_set_vec==ii, 2) - act_frac_sys(frac_set_vec==ii, 4)).^2).^(1/2) );
end

comp_segm_list = [];
ids = find(components_network_nodes(:, max_id) == max_id);
comp_segm_mat = incidence_matrix(ids, :);
for jj = 1:size(comp_segm_mat, 1)
    comp_segm_list = [comp_segm_list, find(comp_segm_mat(jj, :))];
end
comp_segm_list = unique(comp_segm_list);

sets_in_cluster = frac_set_vec(comp_segm_list);
unq_sets_clusters = unique(sets_in_cluster);
length_frac_sets_in_cluster = length_frac_sets(unq_sets_clusters);

[sorted_length, ids_len_set] = sort(length_frac_sets_in_cluster, 'descend');

total_length_cluster = sum(sorted_length);
current_ratio = 1;
iter = 1;
total_remove_segm = [];

max_iter = 4;
if length(unq_sets_clusters) > 55
    max_iter = 33; 
elseif length(unq_sets_clusters) > 45
    max_iter = 28; 
elseif length(unq_sets_clusters) > 35
    max_iter = 23; 
elseif length(unq_sets_clusters) > 25
    max_iter = 17; 
elseif length(unq_sets_clusters) > 15
    max_iter = 8; 
elseif length(unq_sets_clusters) > 10
    max_iter = 6;
end

while (current_ratio - threshold_conn > reduction_in_length) && (iter < length(unq_sets_clusters)) && (iter < max_iter) && (size(new_segm_to_add, 1) == 0)
    % Set 1:
    set_1_comb = unq_sets_clusters(ids_len_set(1:(end - iter)));

    % Set 2:
    set_2_comb = unq_sets_clusters(ids_len_set(end - iter + 1));

    % Create list of segments in each set:
    segm_list_set_1_ids = [];
    for ii = 1:length(set_1_comb)
        segm_list_set_1_ids = [segm_list_set_1_ids; comp_segm_list(sets_in_cluster==set_1_comb(ii))'];
    end

    segm_list_set_2_ids = comp_segm_list(sets_in_cluster==set_2_comb)';

    % Check if segments occur in both sets:
    common_segms = [];
    for ii = segm_list_set_2_ids'
        if any(segm_list_set_1_ids == ii)
            common_segms = [common_segms; ii];
        end   
    end

    % Create list of nodes in each set:
    nodes_id_in_set_1 = [];
    for ii = segm_list_set_1_ids'
        node_ids = find(incidence_matrix(:, ii));
        nodes_id_in_set_1 = [nodes_id_in_set_1; node_ids]; 
    end
    nodes_id_in_set_1 = unique(nodes_id_in_set_1);

    nodes_id_in_set_2 = [];
    for ii = segm_list_set_2_ids'
        node_ids = find(incidence_matrix(:, ii));
        nodes_id_in_set_2 = [nodes_id_in_set_2; node_ids]; 
    end
    nodes_id_in_set_2 = unique(nodes_id_in_set_2);

    % Find common nodes among two sets:
    common_nodes = [];
    for ii = nodes_id_in_set_2'
        if any(nodes_id_in_set_1 == ii)
            common_nodes = [common_nodes; ii];
        end   
    end

    % Find segments in set that contain this node:
    segment_to_remove = [];
    for ii = common_nodes'
        ids_segms = find(incidence_matrix(ii, :));

        remove = false;
        for jj = ids_segms
            if any(jj == segm_list_set_2_ids)
                segment_to_remove = [segment_to_remove; jj];
            end
        end
    end

    % Plot results (set 1, set2 without removed, and remove segments)
    set_1_list = segm_list_set_1_ids;
    set_2_list = segm_list_set_2_ids;
    indexes = true(length(segm_list_set_2_ids), 1);
    for ii = segment_to_remove'
        if any(set_2_list==ii)
            indexes(find(set_2_list==ii)) = false;
        end
    end
    seg_2_rem_list = segment_to_remove;
   
    current_ratio = sum(sorted_length(1:(end - iter))) / total_length_cluster;
    iter = iter + 1; 
    total_remove_segm = [total_remove_segm; segment_to_remove];
end

if plot_results && (iter > 1)
    figure();
    plot(act_frac_sys(set_1_list, [1, 3])', act_frac_sys(set_1_list, [2, 4])', 'color', [0, 0, 1])
    hold on
    plot(act_frac_sys(set_2_list(indexes), [1, 3])', act_frac_sys(set_2_list(indexes), [2, 4])', 'color', [0, 0, 0])
    plot(act_frac_sys(seg_2_rem_list, [1, 3])', act_frac_sys(seg_2_rem_list, [2, 4])', 'color', [1, 0, 0])
    hold off
end

if iter > 1
    % Remove segments:
    indexes_keep = true(size(act_frac_sys, 1), 1);
    for ii = total_remove_segm'
        indexes_keep(ii) = false;
    end

    act_frac_sys = act_frac_sys(indexes_keep, :);
    frac_set_vec = frac_set_vec(indexes_keep);
end

% Update act_frac_sys with new segments:
act_frac_sys = round(act_frac_sys * 10^decimals) * 10^(-decimals);

% Make sure segments are sorted from x_min to x_max:
act_frac_sys_dummy = act_frac_sys;
indices = act_frac_sys_dummy(:, 1) > act_frac_sys_dummy(:, 3);
act_frac_sys(indices, 1:2) = act_frac_sys_dummy(indices, 3:4);
act_frac_sys(indices, 3:4) = act_frac_sys_dummy(indices, 1:2);

[act_frac_sys, id_frac_set, ~] = unique(act_frac_sys, 'rows');
frac_set_vec = frac_set_vec(id_frac_set);

if size(new_segm_to_add, 1) > 0
    % Only if segments are added (i.e. connectivity was increased):
    new_segm_to_add = round(new_segm_to_add * 10^decimals) * 10^(-decimals);

    % Make sure segments are sorted from x_min to x_max:
    new_segm_dummy = new_segm_to_add;
    indices = new_segm_dummy(:, 1) > new_segm_dummy(:, 3);
    new_segm_to_add(indices, 1:2) = new_segm_dummy(indices, 3:4);
    new_segm_to_add(indices, 3:4) = new_segm_dummy(indices, 1:2);

    [new_segm_to_add, id_frac_set, ~] = unique(new_segm_to_add, 'rows');

    for ii = 1:size(new_segm_to_add, 1)
        % New coordinate:
        new_coord = [mean(new_segm_to_add(ii, [1, 3])), mean(new_segm_to_add(ii, [2, 4]))];

        % Find which nodes to move:
        ids_left = find( ( (new_segm_to_add(ii, 1) - act_frac_sys(:, 1)).^2 + ...
                           (new_segm_to_add(ii, 2) - act_frac_sys(:, 2)).^2 ).^(1/2) < tolerance_zero);

        % Update nodes:
        if ~isempty(ids_left)
            act_frac_sys(ids_left, [1, 2]) = ones(size(ids_left, 1), 1) * new_coord;
        end

        % Find which nodes to move:
        ids_left = find( ( (new_segm_to_add(ii, 1) - act_frac_sys(:, 3)).^2 + ...
                           (new_segm_to_add(ii, 2) - act_frac_sys(:, 4)).^2 ).^(1/2) < tolerance_zero);

        % Update nodes:
        if ~isempty(ids_left)
            act_frac_sys(ids_left, [3, 4]) = ones(size(ids_left, 1), 1) * new_coord;
        end

        % Find
        ids_right = find( ( (new_segm_to_add(ii, 3) - act_frac_sys(:, 3)).^2 + ...
                            (new_segm_to_add(ii, 4) - act_frac_sys(:, 4)).^2 ).^(1/2) < tolerance_zero);

        % Update nodes:
        if ~isempty(ids_right)
            act_frac_sys(ids_right, [3, 4]) = ones(size(ids_right, 1), 1) * new_coord;
        end

        % Find
        ids_right = find( ( (new_segm_to_add(ii, 3) - act_frac_sys(:, 1)).^2 + ...
                            (new_segm_to_add(ii, 4) - act_frac_sys(:, 2)).^2 ).^(1/2) < tolerance_zero);

        % Update nodes:
        if ~isempty(ids_right)
            act_frac_sys(ids_right, [1, 2]) = ones(size(ids_right, 1), 1) * new_coord;
        end
    end

    % Update new graph (based on connectivity:
    new_num_comp = num_comp - size(comp_added, 1) + 1;
    new_components_network_segms = zeros(num_main_segm, new_num_comp);

    comp_count = 2;

    for ii = 1:num_comp
        % Find segments belonging to subgraph:
        ids = find(components_network_segms(:, ii) == ii);

        if any(ii == comp_added)
            new_components_network_segms(ids, 1) = 1;
        else
            % Component not in main network:
            new_components_network_segms(ids, comp_count) = comp_count;
            comp_count = comp_count + 1;
        end
    end
    
    if plot_results
        % Plot components of graph:
        colorjet = jet;
        coll_def_size = colorjet(floor(linspace(1, 64, new_num_comp)),:);

        figure()
        hold on

        for ii = 1:new_num_comp
            ids = find(new_components_network_segms(:, ii) == ii);
            plot(act_frac_sys(ids, [1, 3])', act_frac_sys(ids, [2, 4])', 'color', coll_def_size(ii, :))    
        end
        hold off
    end

    % Find again intersections in case something got messed up and extract
    % unique and remove zero segments:
    act_frac_sys = round(act_frac_sys * 10^decimals) * 10^(-decimals);

    % Make sure segments are sorted from x_min to x_max:
    act_frac_sys_dummy = act_frac_sys;
    indices = act_frac_sys_dummy(:, 1) > act_frac_sys_dummy(:, 3);
    act_frac_sys(indices, 1:2) = act_frac_sys_dummy(indices, 3:4);
    act_frac_sys(indices, 3:4) = act_frac_sys_dummy(indices, 1:2);

    [act_frac_sys, id_frac_set, ~] = unique(act_frac_sys, 'rows');
    frac_set_vec = frac_set_vec(id_frac_set);

    calc_intersections;

end

act_frac_sys = round(act_frac_sys * 10^decimals) * 10^(-decimals);

% Make sure segments are sorted from x_min to x_max:
act_frac_sys_dummy = act_frac_sys;
indices = act_frac_sys_dummy(:, 1) > act_frac_sys_dummy(:, 3);
act_frac_sys(indices, 1:2) = act_frac_sys_dummy(indices, 3:4);
act_frac_sys(indices, 3:4) = act_frac_sys_dummy(indices, 1:2);

[act_frac_sys, id_frac_set, ~] = unique(act_frac_sys, 'rows');
frac_set_vec = frac_set_vec(id_frac_set);

new_length_segm = ( (act_frac_sys(:, 1) - act_frac_sys(:, 3)).^2 + ...
                    (act_frac_sys(:, 2) - act_frac_sys(:, 4)).^2 ).^(1/2);
                
act_frac_sys(new_length_segm < tolerance_zero, :) = [];    
[~, ~, current_connectivity] = ...
        calc_connectivity_graph(act_frac_sys, plot_results);
end
end