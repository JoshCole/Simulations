fsArray = [2, 4, 8];
alpha = 0.22;

for i = 1:length(fsArray)
    fs = fsArray(i);
    [ h ] = RootRaiseCosine( alpha, fs );
    figure
    plot(-fs:1/fs:fs,h)
end
%%
fsArray = [2, 4, 8];
alpha = 0.22;

for i = 1:length(fsArray)
    fs = fsArray(i);
    [ h ] = RaisedCosine( alpha, fs );
    figure
    plot(-fs:1/fs:fs,h)
end