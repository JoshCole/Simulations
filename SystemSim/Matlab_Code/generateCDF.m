function generateCDF( data )
colList = hsv(3);
% Downlink
figure
for i = 1:3
    x =data.downlink.thru(:,i);
%     x=data.aggregated.DL_thru(:,i);
    A=min(x); B=max(x);
    [n,x]=hist(x,[A:0.05:B]);
    t = cumsum(n)/sum(n);
    plot(x,t,'col',colList(i,:))
    hold on
    grid on
end
title('Throughput DL')
legend('Sector 1','Sector 2','Sector 3')

figure
for i = 1:3
    x =data.downlink.thru_agg(:,i);
%     x=data.aggregated.DL_thru(:,i);
    A=min(x); B=max(x);
    [n,x]=hist(x,[A:0.05:B]);
    t = cumsum(n)/sum(n);
    plot(x,t,'col',colList(i,:))
    hold on
    grid on
end
title('Throughput DL Agg')
legend('Sector 1','Sector 2','Sector 3')

figure
for i = 1:3
    x =data.downlink.snir(:,i);
%     x=data.aggregated.snir(:,i);
    A=min(x); B=max(x);
    [n,x]=hist(x,[A:0.5:B]);
    t = cumsum(n)/sum(n);
    plot(x,t,'col',colList(i,:))
    hold on
    grid on
end
title('SNIR DL')
legend('Sector 1','Sector 2','Sector 3')

figure
for i = 1:3
    x =data.downlink.interference(:,i);
%     x=data.aggregated.snir(:,i);
    A=min(x); B=max(x);
    [n,x]=hist(x,[A:0.5:B]);
    t = cumsum(n)/sum(n);
    plot(x,t,'col',colList(i,:))
    hold on
    grid on
end
title('Interference DL')
legend('Sector 1','Sector 2','Sector 3')


% Uplink
figure
for i = 1:3
    x =data.uplink.thru(:,i);
%     x=data.aggregated.DL_thru(:,i);
    A=min(x); B=max(x);
    [n,x]=hist(x,[A:0.05:B]);
    t = cumsum(n)/sum(n);
    plot(x,t,'col',colList(i,:))
    hold on
    grid on
end
title('Throughput uL')
legend('Sector 1','Sector 2','Sector 3')

figure
for i = 1:3    
    x =data.uplink.thru_agg(:,i);
%     x=data.aggregated.DL_thru(:,i);
    A=min(x); B=max(x);
    [n,x]=hist(x,[A:0.05:B]);
    t = cumsum(n)/sum(n);
    plot(x,t,'col',colList(i,:))
    hold on
    grid on
end
title('Throughput uL Agg')

legend('Sector 1','Sector 2','Sector 3')

figure
for i = 1:3
    x =data.uplink.snir(:,i);
%     x=data.aggregated.snir(:,i);
    A=min(x); B=max(x);
    [n,x]=hist(x,[A:0.5:B]);
    t = cumsum(n)/sum(n);
    plot(x,t,'col',colList(i,:))
    hold on
    grid on

end
title('SNIR UL')
legend('Sector 1','Sector 2','Sector 3')

figure
for i = 1:3
    x =data.uplink.interference(:,i);
%     x=data.aggregated.snir(:,i);
    A=min(x); B=max(x);
    [n,x]=hist(x,[A:0.5:B]);
    t = cumsum(n)/sum(n);
    plot(x,t,'col',colList(i,:))
    hold on
    grid on
end
title('Interference UL')
legend('Sector 1','Sector 2','Sector 3')

end

