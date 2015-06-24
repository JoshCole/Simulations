close all 
clear all
x = 0.1:0.1:2;          % input power vector
xdB = -15:1:3;
x = 10.^(xdB / 10); 
n = length(x);          % length of x
backoff = 0.0; % backoff
[AM PM y] = salehs_model(x,backoff,n); % nonlinearity model
figure
subplot(2,1,1)
pin = 10*log10(abs(x)); % input power in dB
pout = 10*log10(abs(y)); % output power in dB
plot(pin,pout); grid;
xlabel('Input power - dB')
ylabel('Output power - dB')
subplot(2,1,2)
plot(pin,(180/pi)*unwrap(angle(y))); grid;
xlabel('Input power - dB')

figure
plot(xdB,10*log10(abs(y./x)))
grid on

figure
plot(pin,pin,'g')
hold on
plot(pin,pin-1)
plot(pin,pout,'r')
legend('Theoretical','1 dB Change', 'Amplifier Output')
%%
backoff = input('Enter backoff in dB > ');
f1 = -1.0; f2 = 2.0; ts = 1.0/128; n = 1024;
for k=1:n
t(k) = (k-1)*ts;
x(k) = exp(1i*2*pi*f1*t(k))+0.707*exp(1i*2*pi*f2*t(k));
y(k) = salehs_model(x(k),-1*backoff,1);
end
[psdx,freq] = log_psd(x,n,ts);
[psdy,freq] = log_psd(y,n,ts);
figure
subplot(2,1,1)
plot(freq,psdx); grid; title('Input to the NL');
ylabel('PSD in dB');
subplot(2,1,2)
plot(freq,psdy); grid; title('Output of the NL');
ylabel('PSD in dB'); xlabel('Frequency in Hz');
%%
close all
clear all
backoff = 0;
% backoff = input('Enter backoff in dB > ');
f1 = -1.0; f2 = 2.0; ts = 1.0/128; n = 1024;
for k=1:n
t(k) = (k-1)*ts;
x(k) = exp(1i*2*pi*f1*t(k))+0.707*exp(1i*2*pi*f2*t(k));
[ AM_AM(k), AM_PM(k), y(k) ] = AmplifierModel( x(k),-1*backoff );
end
[psdx,freq] = log_psd(x,n,ts);
[psdy,freq] = log_psd(y,n,ts);
figure
subplot(2,1,1)
plot(freq,psdx); grid; title('Input to the NL');
ylabel('PSD in dB');
subplot(2,1,2)
plot(freq,psdy); grid; title('Output of the NL');
ylabel('PSD in dB'); xlabel('Frequency in Hz');

figure

pin = 10*log10(abs(x)); % input power in dB
pout = 10*log10(abs(y)); % output power in dB
plot(pin,pin,'r'); 
hold on
plot(pin,pin-1,'g')
plot(pin,pout); grid;
xlabel('Input power - dB')
ylabel('Output power - dB')
legend('Input','1 dB Change', 'Output')


figure
plot(pin,10*log10(abs(y./x)))
grid on
%%
N = 512;

M = 4;
[ k, Es, Esnorm, Eb, Ebnorm, SymbArray ] = GetSymbolArrayData( M );
numBits = k*N;

Bits = rand(1,numBits) > 0.5;

% Generate Symbol Index
[ Idx, BlockSize ] = SymbolIndex( Bits, k );
        
% Generate Symbols
Symbols = SymbArray(Idx);

OversampleRate = 2;                                             % Oversample rate
alpha = 0.22;                                                   % Roll Of Rate
truncation = 32;                                                % Truncation in respect to number of symbols
[ InterperlationFilter ] = RaisedCosine( alpha, truncation ,OversampleRate );      % Interperlation Filter

% Oversample
OVA = zeros(1,OversampleRate*N);
OVA(1:OversampleRate:end) = Symbols;

% Apply Interperation Filter
x = conv(OVA,InterperlationFilter);

backoff = input('Enter backoff in dB > ');
ts = 1/OversampleRate;
n = length(x);

for k=1:n
[ AM_AM(k), AM_PM(k), y(k) ] = AmplifierModel( x(k),-1*backoff );
end
[psdx,freq] = log_psd(x,n,ts);
[psdy,freq] = log_psd(y,n,ts);
figure
subplot(2,1,1)
plot(freq,psdx); grid; title('Input to the NL');
ylabel('PSD in dB');
subplot(2,1,2)
plot(freq,psdy); grid; title('Output of the NL');
ylabel('PSD in dB'); xlabel('Frequency in Hz');

figure

pin = 10*log10(abs(x)); % input power in dB
pout = 10*log10(abs(y)); % output power in dB
plot(pin,pin,'r'); 
hold on
plot(pin,pin-1,'g')
plot(pin,pout); grid;
xlabel('Input power - dB')
ylabel('Output power - dB')
legend('Input','1 dB Change', 'Output')


figure
plot(pin,10*log10(abs(y./x)))
grid on
%%
