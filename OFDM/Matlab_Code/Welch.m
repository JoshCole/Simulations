function [ Pxx ] = Welch( x, Nfft, Fs )
M = length(x);
% Hanning Window
% -------------- 
[ win ] = Hanning( Nfft );  
Ls = length(win);

noverlap = fix(0.5.*Ls);

k = (M-noverlap)./(Ls-noverlap);
k = fix(k);

% Define Overlap  
% --------------
LminusOverlap = Ls-noverlap;
xStart = 1:LminusOverlap:k*LminusOverlap;
xEnd   = xStart+Ls-1;

Sxx = zeros(Nfft,1);                    % Power Spectrum Density accumulator
WindowEnergy = sum(abs(win).^2);        % Compensates for the power of the window, DC gain of window function
coherentPowerGain = sum(abs(win)).^2;
[ Normx, ux, Ex, Px, XRMS, simga2, DCGain, PowerGain ] = SignalInfo( win );
for jj = 1:k,
    xk = x(xStart(jj):xEnd(jj));
    xw = xk.*win;
    Xx = fft(xw,Nfft); 
    P = Xx.*conj(Xx)/Px;                % Auto spectrum.
    Sxx = Sxx + P;                      % Sum scaled periodograms
end
Sxx = Sxx./k;                           % Average the sum of the periodogram
Pxx = Sxx./Fs;                          % Scale by the sampling frequency to obtain the psd

MeanPxx = mean(Pxx); 


end

