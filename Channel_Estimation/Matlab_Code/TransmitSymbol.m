function [ RxSymbol, Noise ,No, NoiseVar ] = TransmitSymbol( TxSymbol, Es, Esnorm, EsNo_lin, Fs )
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
FrameSize = length(TxSymbol);
No = (Es./EsNo_lin);    % Noise Spectral Density
NoiseVar = (No/2);      % Noise Variance
AWGN = (((randn(1,FrameSize))+1i*(randn(1,FrameSize)))/(sqrt(2)));     % Normalised AWGN in respect to Esnorm
Noise = (sqrt(No).*AWGN)/sqrt(Fs*Esnorm);                              % Noise to be induced to Signal
RxSymbol = TxSymbol + Noise;    % Recieved Signal
end

