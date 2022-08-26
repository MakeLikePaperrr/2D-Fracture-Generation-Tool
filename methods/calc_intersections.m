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
% Allocate large array for new points:
n_fracs = size(act_frac_sys, 1);
max_frac_int = round( n_fracs / 2 );
max_new_pts = round(max_frac_int * n_fracs * 0.01);
max_length_new_segm = n_fracs + max_new_pts * 2;

% Segment list global:
segment_id = frac_set_vec;
new_segment_id = zeros(max_length_new_segm, 1);

new_points = zeros(max_new_pts, 2);
ith_jj = zeros(max_new_pts, 1);
new_fract_sys = zeros(max_length_new_segm, 4);

ith_pt = 0;
glob_segm_count = 0;

for ii = 1:n_fracs
    % Pointer to first newly added point (if added!)
    ith_old = ith_pt + 1;
    
    for jj = (ii + 1):n_fracs
        ii_frac = act_frac_sys(ii, :);
        jj_frac = act_frac_sys(jj, :);
        
        [t, s, int_coord] = find_parametric_intersect(ii_frac, jj_frac);
        
        if (t >= (0 - tolerance_intersect) && t <= (1 + tolerance_intersect)) && ...
           (s >= (0 - tolerance_intersect) && s <= (1 + tolerance_intersect))
            % Only store intersections of segments that don't already share
            % node: 
            if ~(norm(ii_frac([1, 2]) - jj_frac([1, 2])) < tolerance_intersect || ...
                 norm(ii_frac([1, 2]) - jj_frac([3, 4])) < tolerance_intersect || ...
                 norm(ii_frac([3, 4]) - jj_frac([1, 2])) < tolerance_intersect || ...
                 norm(ii_frac([3, 4]) - jj_frac([3, 4])) < tolerance_intersect)
                % ii and jj intersect:
                ith_pt = ith_pt + 1;

                % Store intersection coordinate:
                new_points(ith_pt, :) = int_coord;

                % Store jj intersection:
                ith_jj(ith_pt) = jj;
            end
        end
    end
    
    % Append all newly found intersections to full array:
    % four extra segments: start, end, new_ii_int, prev_jj_int
    new_ii_int = new_points(ith_old:ith_pt, :);
    prev_jj_int = new_points(ith_jj==ii, :);

    tot_new_pts = 2 + size(new_ii_int, 1) + size(prev_jj_int, 1);
    
    % Start:
    tot_loc_pts_list = act_frac_sys(ii, [1, 2]);
    
    % End: 
    tot_loc_pts_list = [tot_loc_pts_list; act_frac_sys(ii, [3, 4])];
    
    % new_ii_int:
    tot_loc_pts_list = [tot_loc_pts_list; new_ii_int];
    
    % new_ii_int:
    tot_loc_pts_list = [tot_loc_pts_list; prev_jj_int];
    
    % Sort local points from x_min to x_max:
    [unique_x, ids] = unique(tot_loc_pts_list(:, 1));
    
    if length(ids) ~= length(tot_loc_pts_list(:, 1))
        % All x-coordinates are the same to some precision, sort y
        % dimension instead:
        tot_loc_pts_list = sortrows(tot_loc_pts_list,2);
        
    else
        % All are unique so sorting results in actual sort:
        tot_loc_pts_list = sortrows(tot_loc_pts_list,1);
    end
    
    % Loop over points to create segments:
    % kk = 1:(tot_new_pts - 1)
    %   kk --> segm_kk = [x_k, y_k, x_k+1, y_k+1]
    tot_new_segm = tot_new_pts - 1;
    tot_loc_segm_list = zeros(tot_new_segm, 4);
    
    for kk = 1:tot_new_segm
        tot_loc_segm_list(kk, :) = [tot_loc_pts_list(kk, 1), ...
                                    tot_loc_pts_list(kk, 2), ...
                                    tot_loc_pts_list(kk+1, 1), ...
                                    tot_loc_pts_list(kk+1, 2)];
    end
    
    % Store segment list in global array:
    new_fract_sys((glob_segm_count + 1):(glob_segm_count + tot_new_segm), :) = tot_loc_segm_list;
    
    % Store new segment ids:
    new_segment_id((glob_segm_count + 1):(glob_segm_count + tot_new_segm)) = segment_id(ii);
    
    % Update global segment count:
    glob_segm_count = glob_segm_count + tot_new_segm;
end

% Extract only actual segments
act_frac_sys = new_fract_sys(1:glob_segm_count, :);
new_segment_id = new_segment_id(1:glob_segm_count);

% Determine length of new "main" segments:
len_segm_new = sqrt( (act_frac_sys(:, 1) - act_frac_sys(:, 3)).^2 + ...
                     (act_frac_sys(:, 2) - act_frac_sys(:, 4)).^2 );

% Remove non-zero fracture segments:                      
nonzero_segm = len_segm_new > tolerance_zero;
act_frac_sys = act_frac_sys(nonzero_segm, :); 

new_segment_id = new_segment_id(nonzero_segm);
frac_set_vec = new_segment_id;

