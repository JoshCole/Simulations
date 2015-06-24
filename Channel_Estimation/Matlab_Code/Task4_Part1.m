close all
clear all

% M Sequence
G = [1 0 0 0 0 1 1];
GF = G(2:end); 
N = length(G)-1; 
M = 2^N-1;
Reg = zeros(1,N); 
x = zeros(M,1);
% Number of zeros

% Initalise Reg 
Reg(end)=1;   
 
for n = 1:M
    buffer = mod(sum(GF.*Reg),2);      % Caluculate New Bit
    Reg = [buffer,Reg(1:N-1)];         % Shift Register   
    x(n) = Reg(N);
end
% numZeros = find(x == 0);
% numOnes = find(x == 1);
% numZeros = length(numZeros);
% numOnes = length(numOnes);
% if (numOnes-numZeros)== 1
%     disp('Balanced')
% else
%     disp('Unbalanced')
% end

% Correlation

xNew = x.*2-1;
x = xNew;
h = x;
Nx = length(x);
Nh = length(h);
  
N = Nx;       % For Circular Correlation Nx must be equal to Nh
% Time Domain Circular Correlation
y = zeros(1,N);
for n=1:N 
    y(n) = sum(x(mod(-n+1:1:Nh-n,Nx)+1).*h(1:Nh)); 
    % mod(Nk,Nh) if Nh ~= 0,returns Nk + n.*Nh where n = floor(Nk./Nh).
end

figure
stem(abs(y))
title('Circular Correlation')
xlabel('Time');
ylabel('Magnitude');
axis([0 length(y) min(abs(y)) max(abs(y))])


% Frequency Domain Correlation
X = fft(xNew);
AutoCol = X.*conj(X);%abs(X).^2;

% Time Domain Correlation
TimeDomainConv = ifft(AutoCol);
subplot(2,1,2)
stem(TimeDomainConv)
title('Caclulated in Time Domain')

% Time Domain Linear Correlation
for n=1:N  
    xk= mod(-Nh-n:1:-Nh-1,Nx)+1;
    y2(n) = sum(x(xk).*h(1:length(xk))); 
end

% y=zeros(1,N);
% X = x; 
% H = h;
% for n=-Nx:1:Nx;
%     n2=mod(n,Nx)+1; 
%     for k=1:1:n2;
%         y(n2)=y(n2)+X(n+k+1)*H(k);
%     end;
% end;

figure
stem(y2)
h = conj(flipdim([h;h],1));
y = conv(x,h);
figure
plot(abs(y));
title('Linear Correlation')
xlabel('Time');
ylabel('Magnitude');
axis([0 length(y) min(abs(y)) max(abs(y))])

%%
% Convolution
x =[0 0 0 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0] ;
h = [ 0 15 15 15 15 10 10 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0] ;
Nx = length(x); 
Nh = length(h); 
  
N = Nx;       % For Circular convolution Nx must be equal to Nh

% Time Domain Circular Convolution
y = zeros(1,n);
for n=1:Nx 
    y(n) = sum(x(1:Nx).*h(mod(n-1:-1:n-Nh, Nh)+1)); 
    % mod(Nk,Nh) if Nh ~= 0,returns Nk - n.*Nh where n = floor(Nk./Nh).
end
figure
plot(y)

% Time Domain Linear Convolution
% h = [ 0 15 15 15 15 10 10 10 0 0 0 0 0 0 0 ] ;
Nx = length(x); 
Nh = length(h); 

figure
subplot(2,1,1)
y = conv(x,h);
plot(y)

N = Nx+Nh-1;
% Zero Padding of signals x and h
X=zeros(1,N);
H=zeros(1,N);
X(1:Nx) = x;
H(1:Nh) = h;
y=zeros(1,N);

for n=1:1:N;
    for k=1:1:n;
        y(n)=y(n)+X(k)*H(n-k+1);
    end;
end;
subplot(2,1,2)
plot(y)

