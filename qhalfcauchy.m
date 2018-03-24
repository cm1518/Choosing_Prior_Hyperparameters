function x=qhalfcauchy(p,scale)
% evaluate half-cauchy cdf with scale parameter 
% translated from R toolbox "LaplacesDemon"
if(any(p < 0) || any(p > 1)) 
    error('p must be in [0,1].')
end;

if(any(scale <= 0)) 
    error('The scale parameter must be positive.')
end

x = scale*tan((pi*p)/2);
 