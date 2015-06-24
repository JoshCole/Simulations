function [ UpSample, Nsample ] = OverSample( Signal, NyquistRate )
Nsig = length(Signal);
N = NyquistRate*Nsig;
UpSample = zeros(1,N);
UpSample(1:NyquistRate:N) = Signal;
Nsample = length(UpSample);


end

