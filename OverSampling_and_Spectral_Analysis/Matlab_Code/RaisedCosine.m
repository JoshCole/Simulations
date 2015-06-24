function [ h ] = RaisedCosine( alpha, taps ,fs )
t = taps/2;
t_Tc =-t:1/fs:t;

%Raised Cosine Filter

L=length(t_Tc); %Filter Length
R=1; %Data Rate = 1Mbps
T=1/R;
Ts=1/fs;

%Raised Cosing Filter Design

%----------------------------------------------------------

if mod(L,2)==0
    M=L/2 ; % for even value of L
else
    M=(L-1)/2; % for odd value of L
end

h=zeros(1,L); %Place holder for RC filter's transfer function

for n=-M:M
    num=sin(pi*n*Ts/T)*cos(alpha*pi*n*Ts/T);
    den=(pi*n*Ts/T)*(1-(2*alpha*n*Ts/T)^2);

    h(n+M+1)=num/den;        

    if (1-(2*alpha*n*Ts/T)^2)==0
         h(n+M+1)=pi/4*sin(pi*n*Ts/T)/(pi*n*Ts/T);
    end
    if  n==0
         h(n+M+1)=cos(alpha*pi*n*Ts/T)/(1-(2*alpha*n*Ts/T)^2);
    end
end
%-----
end
