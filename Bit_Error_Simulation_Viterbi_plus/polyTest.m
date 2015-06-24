genPoly1 = [1 0 1 ];                 % Polynomials for 3GPP
genPoly2 = [1 1 1];
genPoly = [genPoly1;genPoly2];
nPacket = [1 0 1 1 1 0 0 0];
[n, m] = size(genPoly);            % Number and Length of Polynomials
m= m-1;
k= 1;
GF= zeros(m-1,n);
for i = 1:n
    GF(:,i) = genPoly(i,2:end-1);
end
    reg = zeros(1,m);
    
for i = 1:length(nPacket)
    if reg(1) == 0
        reg = [reg(2:n),nPacket(i)];
    else
        for ii = 2:n
            reg(1:n-1)=xor(GF(ii),reg(2:n));
            reg(n)= xor(1, nPacket(i));
            genOut(i,ii) = reg(n);
        end
    end
end