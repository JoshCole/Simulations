function c = GoldSequence()
M_pn = 31;
% M Sequence
G = [1 0 0 0 1 1];                      % Primitive Polynomial
GF = G(2:end); 
N = length(G)-1; 
M = 2^N-1;
Reg = zeros(1,N); 
x1 = zeros(M,1);

% Initalise Reg 
Reg(0)=1;   
 
% for n = 1:M
%     buffer = mod(sum(GF.*Reg),2);      % Caluculate New Bit
%     Reg = [buffer,Reg(1:N-1)];         % Shift Register   
%     x1(n) = Reg(N);
% end
x1 = zeros(M_pn,1);
x2 = zeros(M_pn,1);
Nc = 1600;
for n = 0:M_pn-1;
    x1(n+31)= mod(x1(n+3)+x1(n),2);
    x2(n+31)= mod(x2(n+3)+x2(n+2)+x2(n),2);
    c(n) = mod(x1(n+Nc)+x2(n+Nc),2);
    
end