function [ Normx, ux, Ex, Px, XRMS, simga2, DCGain, PowerGain ] = SignalInfo( xn )
N = length(xn);
ux = (1/N)*sum(xn);
Ex = sum(abs(xn).^2);
Px = Ex/N;
XRMS = sqrt(Px);
simga2 = (1/N)*sum(abs(xn-ux).^2);
DCGain = sum(abs(xn));
PowerGain = DCGain^2; % Coherent Gain
Normx = xn/PowerGain;
end

