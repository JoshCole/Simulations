% Create input constellation
backoff = input('Enter backoff in dB > ');
N = 1024; % number of points
x1 = 2*fix(4*rand(1,N))-3; % direct components
x2 = 2*fix(4*rand(1,N))-3; % quadrature components
y = x1+1i*x2; % signal space points
%
% Run it thru Saleh’s model
% z = salehs_model(y,-1*backoff,1024);
[ AM_AM, AM_PM, z ] = AmplifierModel( y,-1*backoff );
subplot(1,2,1)
plot(real(y),imag(y));grid; title('Input Constellation');
xlabel('direct'); ylabel('quadrature')
axis equal
subplot(1,2,2)
plot(real(z),imag(z));grid; title('Output Constellation');
xlabel('direct'); ylabel('quadrature')
axis equal
% End of script file.