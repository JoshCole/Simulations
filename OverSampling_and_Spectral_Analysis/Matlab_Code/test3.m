close all
clear all
Fs = 1000;   t = 0:1/Fs:.296;
x = cos(2*pi*t*200)+randn(size(t));  % A cosine of 200Hz plus noise
pwelch(x,[],[],[],Fs,'twosided'); % Uses default window, overlap NFFT. 