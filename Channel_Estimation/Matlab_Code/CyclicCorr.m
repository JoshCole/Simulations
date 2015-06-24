function [ yFreq, yTime ] = CyclicCorr( x,h )
Nx =length(x);
Nh =length(h);

% Add Zero Padding if Required
if (Nx>Nh)
    h1 = [h zeros(1,Nx-Nh)];
else
    h1 = h;
end
if(Nh>Nx)
    x1 = [x zeros(1,Nh-Nx)];
else
    x1 = x;
end


% Frequency Domain

X = fft(x1);
H = fft(h1);

% Cyclic Correlation
yFreq = X.*conj(H);
yTime = ifft(yFreq);

end

