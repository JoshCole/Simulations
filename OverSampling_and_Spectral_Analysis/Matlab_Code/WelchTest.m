clear all
close all
M = 32;
Ks=[1 8 32 128];
nkases = length(Ks);
for kase = 1:4
  K = Ks(kase);
  N = M*K;
  Nfft = 2*M; % zero pad for acyclic autocorrelation
  Sv = zeros(Nfft,1); % PSD "accumulator"
  for m=1:K
    v = randn(M,1);  % noise sample
    V = fft(v,Nfft);
    Vms = abs(V).^2; % same as conj(V) .* V
    Sv = Sv + Vms;   % sum scaled periodograms
  end

  Sv = Sv/K;     % average of all scaled periodograms
  rv = ifft(Sv); % Average Bartlett-windowed sample autocor.
  rvup = [rv(Nfft-M+1:Nfft)',rv(1:M)']; % unpack FFT
  rvup = rvup/M; % Normalize for no bias at lag 0
  figure
  stem(Sv)
  pause
    figure
  stem(rv)
  pause
    figure
  stem(rvup)
  pause
end