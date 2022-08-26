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
% Make sure segments are sorted from x_min to x_max:
act_frac_sys_dummy = act_frac_sys;
indices = act_frac_sys_dummy(:, 1) > act_frac_sys_dummy(:, 3);
act_frac_sys(indices, 1:2) = act_frac_sys_dummy(indices, 3:4);
act_frac_sys(indices, 3:4) = act_frac_sys_dummy(indices, 1:2);

[act_frac_sys, id_frac_set, ~] = unique(act_frac_sys, 'rows');
frac_set_vec = frac_set_vec(id_frac_set);

% Determine length of new "main" segments:
len_segm_new = sqrt( (act_frac_sys(:, 1) - act_frac_sys(:, 3)).^2 + ...
                     (act_frac_sys(:, 2) - act_frac_sys(:, 4)).^2 );

% Remove non-zero fracture segments:                      
nonzero_segm = len_segm_new>tolerance_zero;
act_frac_sys = act_frac_sys(nonzero_segm, :); 
frac_set_vec = frac_set_vec(nonzero_segm);

num_main_segm = size(act_frac_sys, 1);
removed_segm = zeros(num_main_segm, 1);

% 1) Main loop over all segments:
for ith_segm = 1:num_main_segm
    % 2) Check if segment is not removed:
    if ~any(ith_segm == removed_segm)
        % 2.1) Check if segment doesn't have zero length:
        segm_length = sqrt( (act_frac_sys(ith_segm, 1) - act_frac_sys(ith_segm, 3))^2 + ...
                            (act_frac_sys(ith_segm, 2) - act_frac_sys(ith_segm, 4))^2 ); 
        if segm_length == 0
            % Store removed segment:
            removed_segm(ith_segm) = ith_segm;
            act_frac_sys(ith_segm, :) = -999;
            continue;
        end
        % 3) Check for both nodes:
        %---------------------
        %       Left node:
        %---------------------
        % a) Is the node part of another segments?
        coord_left = act_frac_sys(ith_segm, [1, 2]);
        left_int_left = sum(act_frac_sys(:, [1, 2]) == ones(size(act_frac_sys, 1), 1) * coord_left, 2) == 2;
        left_int_right = sum(act_frac_sys(:, [3, 4]) == ones(size(act_frac_sys, 1), 1) * coord_left, 2) == 2;
        
        left_count_left = sum(left_int_left(1:end ~=ith_segm));
        left_count_right = sum(left_int_right(1:end ~=ith_segm));
        left_count_tot = left_count_left + left_count_right;
        
        % b) How many segments share this node?
        if left_count_tot == 1
            % c) What is the angle between the two segments? 
            if left_count_left == 1
                % left node of ith_segm occurs in left row of act_frac_sys
                other_segm = find(left_int_left);
                other_segm = other_segm(other_segm~=ith_segm);
                
            elseif left_count_right == 1
                % left node of ith_segm occurs in right row of act_frac_sys
                other_segm = find(left_int_right);
            end
            vec_ith_segm = [act_frac_sys(ith_segm, 1) - act_frac_sys(ith_segm, 3);
                            act_frac_sys(ith_segm, 2) - act_frac_sys(ith_segm, 4)];
            vec_other_segm = [act_frac_sys(other_segm, 1) - act_frac_sys(other_segm, 3);
                              act_frac_sys(other_segm, 2) - act_frac_sys(other_segm, 4)];
            
            % Calculate actual angle:
            try 
                angle_deg = acos((vec_ith_segm' * vec_other_segm) ...
                            / (norm(vec_ith_segm)*norm(vec_other_segm)) ) * 180 / pi;
            catch 
                a = 0;
            end
     
            % Check angle w.r.t. tolerance:
            if abs(angle_deg) < tol_angle
                % Always check if previous no segments were merged such that 
                % segments overlap (duplicate, non-unique), reason for this
                % might be missed intersections in the discretization process!
                temp_array = act_frac_sys([ith_segm, other_segm], :);
                dummy_array = temp_array;
                indices = temp_array(:, 1) > temp_array(:, 3);
                dummy_array(indices, 1:2) = temp_array(indices, 3:4);
                dummy_array(indices, 3:4) = temp_array(indices, 1:2);
                [row, ids] = unique(dummy_array, 'rows');
                if length(ids) == 2
                    % Segments are unique:
                    % d) Merge segments:
                    if left_count_left == 1
                        % left node of ith_segm occurs in left row of act_frac_sys
                        % therefore take (ohter == right) node of other_segm:
                        % ith_Left --> oth_Right
                        act_frac_sys(ith_segm, [1, 2]) = act_frac_sys(other_segm, [3, 4]);

                    elseif left_count_right == 1
                        % left node of ith_segm occurs in right row of act_frac_sys
                        % therefore take (ohter == left) node of other_segm:
                        % ith_Left --> oth_Left
                        act_frac_sys(ith_segm, [1, 2]) = act_frac_sys(other_segm, [1, 2]);
                    end

                    % Perform small check to see if segment is not collapsed:
                    dummy_segm_length = sqrt( (act_frac_sys(ith_segm, 1) - act_frac_sys(ith_segm, 3))^2 + ...
                                              (act_frac_sys(ith_segm, 2) - act_frac_sys(ith_segm, 4))^2 );
                    if dummy_segm_length < 1e-4
                       a = 0; 
                    end

                    % Store removed segment:
                    removed_segm(other_segm) = other_segm;
                    act_frac_sys(other_segm, :) = -999;
                else
                    % Segments are non-unique (on-top of each other)
                    % Store removed segment:
                    removed_segm(other_segm) = other_segm;
                    act_frac_sys(other_segm, :) = -999;
                end
            end
        end
        
        %---------------------
        %       Right node:
        %---------------------
        % a) Is the node part of another segments?
        coord_right = act_frac_sys(ith_segm, [3, 4]);
        right_int_left = sum(act_frac_sys(:, [1, 2]) == ones(size(act_frac_sys, 1), 1) * coord_right, 2) == 2;
        right_int_right = sum(act_frac_sys(:, [3, 4]) == ones(size(act_frac_sys, 1), 1) * coord_right, 2) == 2;
        
        right_count_left = sum(right_int_left(1:end ~=ith_segm));
        right_count_right = sum(right_int_right(1:end ~=ith_segm));
        right_count_tot = right_count_left + right_count_right;
        
        % b) How many segments share this node?
        if right_count_tot == 1
            % c) What is the angle between the two segments? 
            if right_count_left == 1
                % right node of ith_segm occurs in left row of act_frac_sys
                other_segm = find(right_int_left);
                
            elseif right_count_right == 1
                % right node of ith_segm occurs in right row of act_frac_sys
                other_segm = find(right_int_right);
                other_segm = other_segm(other_segm~=ith_segm);
            end
            vec_ith_segm = [act_frac_sys(ith_segm, 1) - act_frac_sys(ith_segm, 3);
                            act_frac_sys(ith_segm, 2) - act_frac_sys(ith_segm, 4)];
            vec_other_segm = [act_frac_sys(other_segm, 1) - act_frac_sys(other_segm, 3);
                              act_frac_sys(other_segm, 2) - act_frac_sys(other_segm, 4)];
            
            % Calculate actual angle:
            try 
                angle_deg = acos((vec_ith_segm' * vec_other_segm) ...
                            / (norm(vec_ith_segm)*norm(vec_other_segm)) ) * 180 / pi;
            catch
                a = 0;
            end
                    
            % Check angle w.r.t. tolerance:
            if abs(angle_deg) < tol_angle
                % Always check if previous no segments were merged such that 
                % segments overlap (duplicate, non-unique), reason for this
                % might be missed intersections in the discretization process!
                temp_array = act_frac_sys([ith_segm, other_segm], :);
                dummy_array = temp_array;
                indices = temp_array(:, 1) > temp_array(:, 3);
                dummy_array(indices, 1:2) = temp_array(indices, 3:4);
                dummy_array(indices, 3:4) = temp_array(indices, 1:2);
                [row, ids] = unique(dummy_array, 'rows');
                if length(ids) == 2
                    % Segments are unique:
                    % d) Merge segments:
                    if right_count_left == 1
                        % right node of ith_segm occurs in left row of act_frac_sys
                        % therefore take (ohter == right) node of other_segm:
                        % ith_Right --> oth_Right
                        act_frac_sys(ith_segm, [3, 4]) = act_frac_sys(other_segm, [3, 4]);

                    elseif right_count_right == 1
                        % right node of ith_segm occurs in right row of act_frac_sys
                        % therefore take (ohter == left) node of other_segm:
                        % ith_Right --> oth_Left
                        act_frac_sys(ith_segm, [3, 4]) = act_frac_sys(other_segm, [1, 2]);
                    end

                    % Store removed segment:
                    removed_segm(other_segm) = other_segm;
                    act_frac_sys(other_segm, :) = -999;
                else
                    % Segments are non-unique (on-top of each other)
                    % Store removed segment:
                    removed_segm(other_segm) = other_segm;
                    act_frac_sys(other_segm, :) = -999;
                end    
            end
        end
        
        % Do other things:
        a = 0;
    end
end

act_frac_sys = act_frac_sys(removed_segm == 0, :);
frac_set_vec = frac_set_vec(removed_segm == 0);
num_tot_new_segm_reduced = size(act_frac_sys, 1);