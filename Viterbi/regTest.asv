genPoly1 = [1 1 1];
genPoly2 = [1 0 1];

genPoly = [genPoly1;genPoly2];

mem = length(genPoly)-1;

dataStream = [1 0 1 1 1 0 0 0];
for i = 1:2
    temp = mod(conv(dataStream,genPoly(i,:)),2);
    c(i,:) = temp(1:length(dataStream))
end
%Initialize code word to all zeros
cd = zeros(1,2*length(dataStream));

%Assemble code word from n outputs
for i = 0:2-1
    cd(1+i:2:end) = c(i+1,:)
end