function [X,Y] = rand_circ(N,x,y,r)
% Generates N random points in a circle.
% RAND_CIRC(N) generates N random points in the unit circle at (0,0).
% RAND_CIRC(N,x,y,r) generates N random points in a circle with radius r 
% and center at (x,y).
if nargin<2
   x = 0;
   y = 0;
   r = 1;
end
Ns = round(1.28*N + 2.5*sqrt(N) + 100); % 4/pi = 1.2732
X = rand(Ns,1)*(2*r) - r;
Y = rand(Ns,1)*(2*r) - r;
I = find(sqrt(X.^2 + Y.^2)<=r);
X = X(I(1:N)) + x;
Y = Y(I(1:N)) + y;