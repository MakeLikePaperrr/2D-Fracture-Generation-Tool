function [t, s, int_coord] = find_parametric_intersect(ii_frac, jj_frac)
% ii_frac = [x_o, y_o, x_1, y_1]
% jj_frac = [x_o, y_o, x_1, y_1]

% Create linear system:
P_0 = [ii_frac(1); 
       ii_frac(2)];
P = [ii_frac(3) - ii_frac(1); 
     ii_frac(4) - ii_frac(2)];  
 
Q_0 = [jj_frac(1); 
       jj_frac(2)];
Q = [jj_frac(3) - jj_frac(1); 
     jj_frac(4) - jj_frac(2)];

A = zeros(2, 2);
A(:, 1) = P;
       
A(:, 2) = -Q;      

if (abs(det(A)) < 1e-5) || (cond(A) < 1e-5)
    % Singular (parallel lines)
    t = -1;
    s = -1;
    int_coord = [];
    
else
    % Solve system and calculate intersection:
    rhs = P_0 - Q_0;
       
    solution = - A \ rhs;
    
    t = solution(1);
    s = solution(2);
    int_coord = P_0 + t * P;
end