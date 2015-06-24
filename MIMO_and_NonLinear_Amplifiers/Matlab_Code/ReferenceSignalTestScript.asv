close all 
clear all

GroupHop = 0;
m = [1 2];                     % 1<= m >= maximum resource blocks in uplink
N_RB_sc = 12;
u = 0;                                    % Group Number 0,1,...,29
g = 0;                                    % Cyclic Shift Index g = 0,1,2,....12
v=0;
ux = [0:29]; 
colList = hsv(length(ux));
figure
for ii = 1:length(ux)
     
    [ r_u_v, q ] = GenReferenceSignal(N_RB_sc, 2, ux(ii), g, GroupHop, 'PUSCH' );
    subplot(4,1,1)
    plot(abs(r_u_v),'-','col',colList(ii,:))
    hold on
    grid on
    set(gca, 'FontSize',14); 
    xlabel('Index');
    ylabel('Mag');
    title('Time Domain')
    subplot(4,1,2)
    plot(abs(fft(r_u_v)/sqrt(length(r_u_v))),'-','col',colList(ii,:))
    hold on
    grid on
    set(gca, 'FontSize',14); 
    xlabel('Index');
    ylabel('Mag');
    title('Frequency Domain')
    
    [ r_u_v, q ] = GenReferenceSignal(N_RB_sc, 1, ux(ii), g, GroupHop, 'PUSCH' );
    subplot(4,1,3)
    plot(abs(r_u_v),'-','col',colList(ii,:))
    hold on
    grid on
    set(gca, 'FontSize',14); 
    xlabel('Index');
    ylabel('Mag');
    title('Time Domain')
    
    subplot(4,1,4)
    plot(abs(fft(r_u_v)/sqrt(length(r_u_v))),'-','col',colList(ii,:))
    hold on
    grid on
    set(gca, 'FontSize',14); 
    xlabel('Index');
    ylabel('Mag');
    title('Frequency Domain')
end

u1 = [0,5,10, 29]; 

mNum = length(m);
colList = hsv(mNum*length(u1));
count = 1;
figure
subplot(3,1,1)
for i = 1:mNum
    [ r_u_v1, q ] = GenReferenceSignal(N_RB_sc, m(i), u, g, GroupHop, 'PUSCH' );
    for ii = 1:length(u1)
        [ r_u_v2, q ] = GenReferenceSignal(N_RB_sc, m(i), u1(ii), g, GroupHop,'PUSCH' );
        plot(0:length(r_u_v1)-1,abs((r_u_v1.*r_u_v2)),'col',colList(count,:))
        Leg{count} =['M_{RS}_{sc}',num2str(m(i)*12), ' u = 0 ','u1 = ', num2str(u1(ii))];
        count = count+1;
    end

end
grid on
set(gca, 'FontSize',14); 
xlabel('Index');
ylabel('Mag');
legend (Leg)

subplot(3,1,2)
hold on
count = 1;
for i = 1:mNum
    [ r_u_v1, q ] = GenReferenceSignal(N_RB_sc, m(i), u, g, GroupHop,'PUSCH' );
    for ii = 1:length(u1)
        [ r_u_v2, q ] = GenReferenceSignal(N_RB_sc, m(i), u1(ii), g, GroupHop,'PUSCH' );
        plot(0:length(r_u_v1)-1,angle((r_u_v1.*r_u_v2)),'col',colList(count,:))
        Leg{count} =['M_{RS}_{sc}',num2str(m(i)*12), ' u = 0 ','u1 = ', num2str(u1(ii))];
        count = count+1;
    end

end
grid on
set(gca, 'FontSize',14); 
xlabel('Index');
ylabel('Phase');
legend (Leg)

subplot(3,1,3)
hold on
count = 1;
for i = 1:mNum
    [ r_u_v1, q ] = GenReferenceSignal(N_RB_sc, m(i), u, g, GroupHop,'PUSCH' );
    for ii = 1:length(u1)
        [ r_u_v2, q ] = GenReferenceSignal(N_RB_sc, m(i), u1(ii), g, GroupHop,'PUSCH' );
        stem(0:length(r_u_v1)-1,abs(ifft((r_u_v1.*conj(r_u_v2)))),'col',colList(count,:))
        Leg{count} =['M_{RS}_{sc}',num2str(m(i)*12), ' u = 0 ','u1 = ', num2str(u1(ii))];
        count = count+1;
    end

end
grid on
set(gca, 'FontSize',14); 
xlabel('Index');
ylabel('Mag');
legend (Leg)


figure
subplot(2,1,1)
hold on
m = [8];                     % 1<= m >= maximum resource blocks in uplink
u = 0;                                    % Group Number 0,1,...,29
g1 = [0:1:11]; 
g = 0;                                    % Cyclic Shift Index g = 0,1,2,....12
v=0;
mNum = length(m);
colList = hsv(mNum*length(g1));
count = 1;
for i = 1:mNum
    [ r_u_v1, q ] = GenReferenceSignal(N_RB_sc, m(i), u, g, GroupHop, 'PUSCH' );
    for ii = 1:length(g1)
        [ r_u_v2, q ] = GenReferenceSignal(N_RB_sc, m(i), u, g1(ii), GroupHop,'PUSCH' );
        stem(0:length(r_u_v1)-1,abs(ifft((r_u_v1.*conj(r_u_v2)))),'col',colList(count,:))
        Leg{count} =['M_{RS}_{sc}',num2str(m(i)*12), ' g = 0 ','g1 = ', num2str(g1(ii))];
        count = count+1;
    end

end
grid on
set(gca, 'FontSize',14); 
title('PUSCH')
xlabel('Index');
ylabel('Mag');
legend (Leg)

subplot(2,1,2)
hold on
u = 0;                                    % Group Number 0,1,...,29
g1 = [0:1:7]; 
N_RB_sc = 12;
g = 0;                                    % Cyclic Shift Index g = 0,1,2,....12
v=0;
mNum = length(m);
colList = hsv(mNum*length(g1));
count = 1;
for i = 1:mNum
    [ r_u_v1, q ] = GenReferenceSignal(N_RB_sc, m(i), u, g, GroupHop, 'SRS' );
    for ii = 1:length(g1)
        [ r_u_v2, q ] = GenReferenceSignal(N_RB_sc, m(i), u, g1(ii), GroupHop,'SRS' );
        stem(0:length(r_u_v1)-1,abs(ifft((r_u_v1.*conj(r_u_v2)))),'col',colList(count,:))
        Leg{count} =['m_{RS}_{sc}',num2str(m(i)*12/2), ' g = 0 ','g1 = ', num2str(g1(ii))];
        count = count+1;
    end

end
grid on
set(gca, 'FontSize',14); 
title('Sounding Reference Signal')
xlabel('Index');
ylabel('Mag');
legend (Leg)


