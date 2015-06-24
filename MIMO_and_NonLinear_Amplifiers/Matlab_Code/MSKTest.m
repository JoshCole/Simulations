N = 5*10^5; % number of bits or symbols
bt = .3;          % 3-dB bandwidth-symbol time
o = 8;            % Oversampling factor
n = 2;            % 2 symbol periods to the filters peak
h = gaussfir(bt,n,o); 

fsHz = 1; % sampling period
T    = 4; % symbol duration    
ct = cos(pi*[-T:N*T-1]/(2*T));
st = sin(pi*[-T:N*T-1]/(2*T));

% MSK Transmitter
    ipBit = rand(1,N)>0.5; % generating 0,1 with equal probability
    ipMod =  2*ipBit - 1; % BPSK modulation 0 -> -1, 1 -> 0
%     ipMod = conv(ipMod,h);
    
    ai = kron(ipMod(1:2:end),ones(1,2*T));  % even bits
    aq = kron(ipMod(2:2:end),ones(1,2*T));  % odd bits

    ai = [ai zeros(1,T)  ]; % padding with zero to make the matrix dimension match 
    aq = [zeros(1,T) aq ];  % adding delay of T for Q-arm
        
    % MSK transmit waveform
    xt = 1/sqrt(T)*[ai.*ct + j*aq.*st];
    
    figure
    plot(real(xt))
    axis([0 100 -1 1])
    