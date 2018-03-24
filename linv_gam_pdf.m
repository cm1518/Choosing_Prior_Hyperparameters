function f = linv_gam_pdf (x, a,b)
% PURPOSE: returns the log pdf at x of the inv_gamma(a,b) distribution, See
% Wikipedia
%---------------------------------------------------
% USAGE: pdf = gamm_pdf(x,mu,nu)
% where: x = a scalar  
%        a = shape
%        b = scale
%---------------------------------------------------
% RETURNS:
%        a vector of pdf at each element of x of the inv_gamma(a,b) distribution      
% --------------------------------------------------



if nargin ~= 3
error('Wrong # of arguments to inv_gamm_pdf');
end;

if any(any(a<=0))
   error('gamm_pdf: parameter a is wrong')
end

if any(any(b<=0))
   error('gamm_pdf: parameter b is wrong')
end


%f = ((nu-2)/2)*log(x) + (-(x*nu)/(2*mu))-cG;
f =a.*log(b)-gammaln(a)+(-a-1).*log(x)-(b./x);
