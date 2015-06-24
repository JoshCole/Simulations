function [ win ] = Hanning( N )
n =1:N;
win = 0.5*(1-cos((2*pi*n)/(N-1))).';
win = win/(N/2);
% win = win/sum(abs(win).^2);
end

