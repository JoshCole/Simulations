function [ codeWord, v, newFrameSize ] = convCode( msg, g )

msg = double(msg);
N = length(g);
v = cell(1,N);
% Generates output corresponding to each generator g0,g1,---gn
% ------------------------------------------------------------
for i = 1:N
    v{i} = mod(conv(msg,g{i}),2);
end
n = length(v{i});
codeWord = false(1,N*n);

% Adds v cells togethere to create code word
% ------------------------------------------
for ii = 0:N-1
    codeWord(ii+1:N:end) = v{ii+1};
end
[M, newFrameSize] = size(codeWord);
end

