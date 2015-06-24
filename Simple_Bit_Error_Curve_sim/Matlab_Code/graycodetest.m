clear all
txBits = [0;0;1;0];

grayCode = false(length(txBits),1);
grayCode(1,1) = txBits(1,1);

cnt = 1;

while cnt<length(txBits)
    a = txBits(cnt,1)
    b = txBits(cnt+1,1);
    grayCode(cnt+1,1) = mod(a + b,2);
    cnt= cnt+1;
end
