function [ Wn ] = Hanning( N )
Wn = 0.5*(1-cos(2*pi*(1:N)'/(N+1)));
end

