Ny = 4;
[ FilterRRC ] = (RootRaiseCosine( 0.22, Ny ));

FilterRRC = fftshift(abs(fft(FilterRRC,1024)));
F = -Ny/2:Ny/(length(FilterRRC)-1):Ny/2;
F = F/(Ny);
plot(F,FilterRRC);