function [ NoiseArray ] = noiseBlockArray(Eb, EbNo_lin ,AWGN )
[M,N] = size(AWGN);
NoiseArray = cell(M,N);
for i = 1:M
    for ii = 1:N
        [m,n] = size(AWGN{M,N});
        NoiseArray{i,ii} = repmat((sqrt(Eb./(EbNo_lin))),1,n).* AWGN{i};
    end
end
end

