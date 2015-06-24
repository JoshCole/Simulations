function [ RxSymbol, AWGN ,No, NoiseVar, rayleighFade ] = RayleighFade( TxSymbol, Es, EsNo_lin, NyquistRate, Method )
% TransmitSymbol
%
% This function calculates the Additive White Gaussian Noise to be induced
% to the current signal. In addition it calcuates the associated noise
% factors
%
% Usage :
%
% [ RxSymbol, Noise ,No, NoiseVar ] = TransmitSymbol( TxSymbol, Es, Esnorm,
% EsNo_lin, FrameSize )
%
% Where         TxSymbol    = Symbols to be transmitted
%
%				Es          = Energy per Symbol, Unitity due to
%                             Normalisation
%
%				Esnorm      = The Normalisation factor for the current
%                             Symbol array, i.e. for 16-QAM Esnorm = 10
%
%				EsNo_lin	= The signal to noise ratio in linear form
%
% TxSymbol = TxSymbol';
N = length(TxSymbol);
No = (Es./EsNo_lin);    % Noise Spectral Density
NoiseVar = No;      % Noise Variance
AWGN = sqrt(NyquistRate)*sqrt(No)*(randn(1,N)+1i*randn(1,N))*(1/sqrt(2));     % Normalised AWGN in respect to Esnorm
x = randn(1, 1); y = randn(1, 1); % Random Variables X and Y with Gaussian Distribution mean zero variance one
rayleighFade = (x + 1i*y)/sqrt(2);       % Normalised Ralyiegh random Variable

if strcmp(Method,'Rayleigh')
    RxSymbol = rayleighFade*TxSymbol ;    % Recieved Signal
else
    RxSymbol = TxSymbol ;    % Recieved Signal
end

