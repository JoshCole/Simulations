M=[1 0 1 1 0 0 1 0 0 0 0];
G = [1 0 0 1 1];
GF = G(2:end-1);
n=length(G)-1;
reg = zeros(1,n);
for k=1:length(M);
%     if reg(1) == 0
%         reg =circshift(reg,[1 -1]);
%         reg(n) = M(k);
%     else
temp =reg(1);
F=reg(1)*GF;
    reg(1:n-1)=xor(F,reg(2:n));
    reg(n)= xor(temp, M(k));
%     end 
end
