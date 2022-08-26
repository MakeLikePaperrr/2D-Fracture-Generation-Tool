function [ mult_h, fraq_length ] = get_fracture_length_powerlaw( xmin, alpha, h )
fraq_length_help = xmin * (1 - rand).^(-1 / (alpha - 1));
mult_h = round(fraq_length_help / h);
fraq_length= mult_h * h;
end

