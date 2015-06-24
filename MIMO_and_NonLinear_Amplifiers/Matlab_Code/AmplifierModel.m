function [ AM_AM, AM_PM, y ] = AmplifierModel( MUx,backoff )
% PA Characteristics
AlphaA = 1;
BetaA = 1;
AlphaPhi = pi/48;
BetaPhi = 0.5;
p=3;
n = length(MUx);
y = zeros(1,n)*(1.0+1i); % initialize output array
AM_AM = zeros(1,n);
AM_PM = zeros(1,n);
phase = zeros(1,n);

absMUx = abs(MUx);
phaseMUx = angle(MUx);
Bk = 10^(backoff/20);
for k =1:n
    x = Bk*absMUx(k);
    AM_AM(k) = (AlphaA.*x)./(1+(BetaA.*x).^2*p).^(1/(2*p));
    AM_PM(k) = (AlphaPhi.*x.^2)./(1+BetaPhi.*x.^2);  %change in Phase
    phase(k) = phaseMUx(k)+ AM_PM(k);
    y(k) = AM_AM(k).*exp(1i*phase(k));
end



end

