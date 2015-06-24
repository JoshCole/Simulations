function [aout, thetaout y]=salehs_model(x,backoff,n)
% This function implements Saleh’s model
% x is the complex input vector of size n;
% Back-off is in db; the input amplitude is scaled by
% c=10^(backoff/20);
% The maximum normalized input power should be less than 3 dB
% i.e., 20 log10(a*abs(x)) < 3 dB
y = zeros(1,n)*(1.0+1i); % initialize output array
a1=2.1587; b1=1.15; % model parameters
a2=4.0; b2=2.1; % model parameters
c=10^(backoff/20); % backoff in dB
for k=1:n
    ain = c*abs(x(k));
    thetain(k) = angle(x(k));
    aout = a1*ain/(1+b1*ain^2);
    thetapm(k) = a2*ain^4/(1+b2*ain^2);
    thetaout(k) = thetain(k)+thetapm(k);
    y(k) = aout*exp(1i*thetaout(k));
end;
% End of function file.
end