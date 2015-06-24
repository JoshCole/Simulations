A = [ 1 2 3]
B = circshift(A,1) % circularly shifts first dimension values down by 1.
B = circshift(A,[1 -1]) % circularly shifts first dimension values
   
R= delN;
if R~=0&&2*R<=N
    q=N/R;
else
    q=N/(R-N);
end
if q is even
    q'=q+gcd(|q|,F)/F;
else
    q'=q;
end
for x = 0:F-1
    s[x*q',mod(F)]=(||x*q'||div(F))
end