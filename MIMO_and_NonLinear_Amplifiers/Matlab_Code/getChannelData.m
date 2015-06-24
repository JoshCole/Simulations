function [ F1out, F2out, Channel ] = getChannelData( fc , F, data ,n)
if length(fc) ==1
    if fc~=0
        [F1min F1] = min(abs(F-(fc-(2.25))));
    else
        [F1min F1] = min(abs(F-(fc)));
    end
    [F2min F2] = min(abs(F-(fc+(2.25))));
    Channel = data(1,F1:F2-1);
    F1out = F(F1);
    F2out = F(F2);
else
    Channel = [];
    for i = 1:length(fc)
        if fc(i)~=0
            [F1min F1] = min(abs(F-(fc(i)-(2.25))));
        else
            [F1min F1] = min(abs(F-(fc(i))));
        end
        [F2min F2] = min(abs(F-(fc(i)+(2.25))));
        temp = data(1,F1:F2-1);
        F1out(i) = F(F1);
        F2out(i) = F(F2);
        Channel = [Channel temp];
    end
end
    
end

