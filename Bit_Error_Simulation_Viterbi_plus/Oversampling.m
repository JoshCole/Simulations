function [ RxOversample ] = Oversampling( RxSymbol, Rate )
LenRx = length(RxSymbol);
ZeroPadding = zeros(1,Rate-1);
RxOversample = zeros(1, LenRx*Rate);
count_1 = 1;
count_2 = Rate;
for i = 1:LenRx
    
    RxOversample(count_1:count_2) = [RxSymbol(i),ZeroPadding];
    count_1 = count_2+1;
    count_2 = (count_2+Rate);
    
end
end
