function [ minDiff, diffIdx ] = eculideanDistance( symbArray, rxSymb)
[M,N] = size(rxSymb);
minDiff = cell(M,N);
diffIdx = cell(M,N);

for i = 1:M
    for ii = 1:N
        [m,n] = size(rxSymb{i,ii});
        for j = 1:m
            for jj = 1:n
               [min_d, Idx] = min(reshape(abs(symbArray- rxSymb{i,ii}(j,jj)).^2,length(symbArray),[]));
               minDiff{i,ii}(j,jj) = min_d;
               diffIdx{i,ii}(j,jj) = Idx-1;
            end
        end
    end
end

