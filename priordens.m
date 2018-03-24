function lpriordens=priordens(x,options_)

switch options_.prior
    
    case 'inv-gamma' % inverse gamma (scaled-inverse-chi2 specification)
        % mean=scale*df/(df-2) for df >2
        % mode=scale*df/(df+2)
        % variance exist for df >4
        a=options_.df/2;
        b=options_.scale*a;
        lpriordens=linv_gam_pdf (x, a,b);
    case 'half-cauchy'
        % no moments exist
        lpriordens=dhalfcauchy(x,options_.scale);
    case 'half-t'
        lpriordens=dhalft(x,options_.scale,options_.df);
    case 'uniform'
        n=length(x);
        lpriordens=-log(options_.upper_bound-options_.lower_bound)*n;
        
end;