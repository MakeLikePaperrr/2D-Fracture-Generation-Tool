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
function [size_subgraphs, len_subgraphs, conn_ratio] = ...
    calc_connectivity_graph(act_frac_sys, plot_results)

% Some constants:
tolerance_zero = 1e-4;

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
    components_network_nodes(abs(component_basis(:, ii) - 1) < tolerance_zero, :) = ii;
end

for ii = 1:num_comp
    comp_segm_list = [];
    comp_segm_mat = incidence_matrix(components_network_nodes(:, ii) == ii, :);
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

size_subgraphs = component_num_segm;
len_subgraphs = component_segm_tot_length;
conn_ratio = max_connected_ratio;