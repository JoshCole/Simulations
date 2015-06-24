%%
clear all 
close all
backoff = input('Enter backoff in dB > ');
f1 = -1.0; f2 = 2.0; ts = 1.0/128; n = 1024;
GaindB = -15:1:5;
gn = length(GaindB);
count  = 1;
for m = 1:gn
    xnew = zeros(1,n);
    y = zeros(1,n);
    for k=1:n
        
    t(k) = (k-1)*ts;
    G = 10^(GaindB(m)/10);
    x(k) = G.*((exp(1i*2*pi*f1*t(k))+exp(1i*2*pi*f2*t(k)))/sqrt(2));
    [ AM_AM(k), AM_PM(k), y(k) ] = AmplifierModel(x(k),-1*backoff);
    end
    [psdx,freq] = log_psd(x,n,ts);
    [psdy,freq] = log_psd(y,n,ts);
    
%     figure
%     subplot(2,1,1)
%     plot(freq,psdx); grid; title('Input to NL Amplifier');
%     ylabel('PSD in dB');
%     subplot(2,1,2)
%     plot(freq,psdy); grid; title('Output of NL Amplifier');
%     ylabel('PSD in dB'); xlabel('Frequency in Hz');

    idx1 = find(freq==-1);
    idx2 = find(freq==2);
    idx3 = find(freq==-4);
    idx4 = find(freq==5);
    P0(m) = psdy(idx1);
    P03(m) = psdy(idx3);
    IMD3(m) = P0(m)-P03(m);
    OIP3(m) = (IMD3(m)/2)+P0(m);
    

    
    N = length(x);
    Ex = sum(abs(x).^2);
    Px = Ex/N;
    xXRMS(m) = sqrt(Px);
    
    N = length(y);
    Ex = sum(abs(y).^2);
    Px = Ex/N;
    yXRMS(m) = sqrt(Px);
    
    G(m) = 10*log10(yXRMS(m)./xXRMS(m));
    IIP3(m) = OIP3(m) - G(m); 
    grad = diff(P03);
    if~isempty(grad)
    if m ==grad(m-1)<=3;
        temp = P03(m);
        temp2 = xXRMS(m);
        temp3 = yXRMS(m);
        grad1=diff(P03(m-1:m));
        grad2=diff(xXRMS(m-1:m));
        grad3=diff(yXRMS(m-1:m));
        for i = 1:20
            P03extetion(i) = temp;
            inputextetion(i) = temp2;
            outputextetion(i) = temp3;
            temp = temp+3;
            temp2 = temp2+1;
            temp3 = temp3+1;
        end
    end
    end
end



pin = 10*log10(abs(xXRMS)); % input power in dB
pout = 10*log10(abs(yXRMS)); % output power in dB

pinextention = 10*log10(abs(inputextetion)); % input power in dB
poutextention = 10*log10(abs(outputextetion)); % output power in dB

figure
plot(pin,pin)
hold on
plot(pin,pout,'r'); 
plot(pin,P03,'--g'); 
plot(pin,IIP3,'k'); 
xlabel('Input power - dB')
ylabel('Output power - dB')
legend('Input power', 'Output power', 'Output third order product','Intercept Point')
axis([-15 30 -60 40])

% figure
% plot(pin,pout); grid;
% xlabel('Input power - dB')
% ylabel('Output power - dB')
% subplot(2,1,2)
% plot(pin,(180/pi)*unwrap(angle(y))); grid;
% xlabel('Input power - dB')
% 
% figure
% plot(pin,10*log10(abs(y./x)))
% grid on
% xlabel('Input power - dB')
% ylabel('Output power - dB')


%%
close all 
clear all
x = 0.1:0.1:2;          % input power vector
xdB = -15:1:3;
x = 10.^(xdB / 10); 
n = length(x);          % length of x
backoff = 0.0; % backoff
[ AM_AM, AM_PM, y ] = AmplifierModel(x,-1*backoff);
figure
subplot(2,1,1)
pin = 10*log10(abs(x)); % input power in dB
pout = 10*log10(abs(y)); % output power in dB
plot(pin,pout); grid;
xlabel('Input power - dB')
ylabel('Output power - dB')
subplot(2,1,2)
phase = (180/pi)*unwrap(angle(y));
plot(pin,phase); grid;
xlabel('Input power - dB')
ylabel('Phase')

figure
plot(xdB,10*log10(abs(y./x)))
xlabel('Input power - dB')
ylabel('Gain dB')
grid on

figure
plot(pin,pin,'g')
hold on
plot(pin,pin-1)
plot(pin,pout,'r')
xlabel('Input power - dB')
ylabel('Output power - dB')
grid on
legend('Theoretical','1 dB Change', 'Amplifier Output')
%%
backoff = input('Enter backoff in dB > ');
f1 = -1.0; f2 = 2.0; ts = 1.0/128; n = 1024;
for k=1:n
t(k) = (k-1)*ts;
x(k) = (exp(1i*2*pi*f1*t(k))+exp(1i*2*pi*f2*t(k)))/sqrt(2);
[ AM_AM(k), AM_PM(k), y(k) ] = AmplifierModel(x(k),-1*backoff);
end
[psdx,freq] = log_psd(x,n,ts);
[psdy,freq] = log_psd(y,n,ts);
figure
subplot(2,1,1)
plot(freq,psdx); grid; title('Input to NL Amplifier');
ylabel('PSD in dB');
subplot(2,1,2)
plot(psdy); grid; title('Output of NL Amplifier');
ylabel('PSD in dB'); xlabel('Frequency in Hz');

pin = 10*log10(abs(x)); % input power in dB
pout = 10*log10(abs(y)); % output power in dB

figure
plot(pin,pout); grid;
xlabel('Input power - dB')
ylabel('Output power - dB')
subplot(2,1,2)
plot(pin,(180/pi)*unwrap(angle(y))); grid;
xlabel('Input power - dB')

figure
plot(pin,10*log10(abs(y./x)))
grid on
xlabel('Input power - dB')
ylabel('Output power - dB')

idx1 = find(freq==-1);
idx2 = find(freq==2);
idx3 = find(freq==-4);
idx4 = find(freq==5);
P0 = psdy(idx1);
P03 = psdy(idx3);
IMD3 = P0-P03;
OIP3 = (IMD3/2)+P0;
    
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
