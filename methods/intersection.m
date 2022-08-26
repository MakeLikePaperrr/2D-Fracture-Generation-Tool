function out=intersection(x1, y1, x2, y2, x3, y3, x4, y4)
    
m1 = (y2-y1)/(x2-x1);
m2 = (y4-y3)/(x4-x3);

b1 = (-m1*x1)+y1;
b2 = (-m2*x3)+y3;

A = [m1 -1; m2 -1];

% rcond(A)
% det(A)

if (abs(det(A)) < 1e-10) || (abs(rcond(A)) < 1e-10) || ...
        any([m1, m2]== Inf) || isnan(rcond(A))
    out = [];
    
else
    b = [-b1; -b2];

    %u = A^-1 * b;

    x = linsolve(A,b);  %solution of the system of equations two compon (x and y)

    % wikipedia https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection

    if x(1)>min(x1,x2) && x(1)<max(x1,x2) && x(1)>min(x3,x4) && x(1)<max(x3,x4)...
       && x(2)>min(y1,y2) && x(2)<max(y1,y2) && x(2)>min(y3,y4) && x(2)<max(y3,y4)
        out= x;
    else
        out=[];
    end
end