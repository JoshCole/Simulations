function [ Xnorm, Averagepower, AbsolutePower, Power, RMS, Simga2, DCGain, PowerGain ] = SignalInfo( xn )
N = length(xn);
Averagepower = (1/N)*sum(xn);
AbsolutePower = sum(abs(xn).^2);
Power = AbsolutePower/N;
RMS = sqrt(Power);
Simga2 = (1/N)*sum(abs(xn-Averagepower).^2);
DCGain = sum(abs(xn));
PowerGain = DCGain^2; % Coherent Gain
Xnorm = xn/PowerGain;
end

