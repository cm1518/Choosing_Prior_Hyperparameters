function [aols,vaols,ls2] = ts_prior_ar(rawdat,tau,p,const)

yt = rawdat(p+1:tau,:);
ylag = mlag(rawdat,p);
if const==1 % including constant
    ylag = [ones(tau-p,1) ylag(p+1:tau,:)];
else
    ylag = ylag(p+1:tau,:);
end;


aols=(ylag'*ylag)\ylag'*yt;

ehat=yt-ylag*aols;
s2=(ehat'*ehat)/(tau-1-p);
ls2=log(s2);

vaols=s2*inv(ylag'*ylag);