function [xt,yt]=discretize_frac_line(xb,yb,h,mult_h,angle,beg)
    
 if beg == 1 
     xt   = xb + h * cosd(angle) * (0:1:mult_h); % = x2
     yt   = yb + h * sind(angle) * (0:1:mult_h); %yb + fraq_length * sind(angle); % = y2
 else
     rand_int = randi([2 mult_h-1],1); %chosing an integer from the line 1 to mult_h (but from mult_h values: discretized)
     
     xt   = xb + h * cosd(angle) * ((0:1:mult_h)-rand_int); %-rand_int % aqui el punto xb=x_disc1(rand_int). Por tanto xt es la misma formula...
     yt   = yb + h * sind(angle) * ((0:1:mult_h)-rand_int); % pero sustrayendole el rand_int que es el punto de interseccion
 end
 