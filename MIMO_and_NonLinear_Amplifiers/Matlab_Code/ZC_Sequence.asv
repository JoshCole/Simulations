close all
m = [25 29 34];
Nrs_zc = 63;
colList = hsv(length(m));
figure
for i = 1:length(m)
    for n = 0:Nrs_zc-1
        ZadoffChu(n+1) = exp(-1i*(pi*m(i)*n*(n+1))/Nrs_zc);
        Leg{i} = ['Root ', num2str(m(i))];
    end
    h = conj(flipdim(ZadoffChu,1));
    y = conv(ZadoffChu,h);
    subplot(2,1,1)
    plot((fft(ZadoffChu)*1/sqrt(Nrs_zc)),'o','col',colList(i,:))
    grid on
    set(gca, 'FontSize',14); 
    xlabel('I');
    ylabel('Q');
    legend(Leg)
    hold on
    subplot(2,1,2)
    hold on
    [F T]=CyclicCorr(ZadoffChu,ZadoffChu);
    stem(real(T),'col',colList(i,:))
    grid on
    set(gca, 'FontSize',14); 
    xlabel('Index');
    ylabel('Magnitude');
    legend(Leg)
    
end 
figure
Ncp = 4;
shift = 63-Ncp;

for i = 1:length(m)
    for n = 0:Nrs_zc-1
        ZadoffChu(n+1) = exp(-1i*(pi*m(i)*n*(n+1))/Nrs_zc);
        Leg{i} = ['Root ', num2str(m(i))];
    end
    
    subplot(2,1,1)
    [F T]=CyclicCorr(ZadoffChu,ZadoffChu);
    stem(0:Nrs_zc-1,real(T),'col',colList(i,:))
    grid on
    set(gca, 'FontSize',14); 
    xlabel('Index');
    ylabel('Magnitude');
    legend(Leg)
    axis([-1 10 0 63])

    hold on
    
    subplot(2,1,2)
    endpoint = Nrs_zc-shift;
    temp = ZadoffChu(endpoint+1:end);
    temp2 = ZadoffChu(1:endpoint);
    ZC_shift = [temp temp2];
    [F T]=CyclicCorr(ZadoffChu,ZC_shift);
    stem(0:Nrs_zc-1,real(T),'col',colList(i,:))
    grid on
    set(gca, 'FontSize',14); 
    xlabel('Index');
    ylabel('Magnitude');
    legend(Leg)
        axis([-1 10 0 63])
    hold on
    
end 
figure
for i = 1:length(m)
    for n = 0:Nrs_zc-1
        ZadoffChu(n+1) = exp(-1i*(pi*m(i)*n*(n+1))/Nrs_zc);
        Leg{i} = ['Root ', num2str(m(i))];
    end
    h = conj(flipdim(ZadoffChu,1));
    y = conv(ZadoffChu,h);
    subplot(3,1,1)
    plot((abs(fft(ZadoffChu))*1/sqrt(Nrs_zc)),'-','col',colList(i,:))
    grid on
    set(gca, 'FontSize',14); 
    xlabel('I');
    ylabel('Q');
%     axis([-1 1 -1 1])
    title('Frequency Domain')
    legend(Leg)
    hold on
    
    subplot(3,1,2)
    plot(abs(ZadoffChu),'-','col',colList(i,:))
    grid on
    set(gca, 'FontSize',14); 
    xlabel('I');
    ylabel('Q');
%     axis([-1 1 -1 1])
    title('Time Domain')
    legend(Leg)
    hold on
    
    subplot(3,1,3)
    hold on
    [F T]=CyclicCorr(ZadoffChu,ZadoffChu);
    stem(real(T),'col',colList(i,:))
    grid on
    set(gca, 'FontSize',14); 
    xlabel('Index');
    ylabel('Magnitude');
    title('Auto Correlation')
    legend(Leg)
    
end 
% % Cross Correlation
% figure
% Ncp = 0;
% shift = 63-Ncp;
% 
% for i = 1:length(m)
%     for n = 0:Nrs_zc-1
%         ZadoffChu_1(n+1) = exp(-1i*(pi*25*n*(n+1))/Nrs_zc);
%         ZadoffChu(n+1) = exp(-1i*(pi*m(i)*n*(n+1))/Nrs_zc);
%         Leg{i} = ['Root ', num2str(m(i))];
%     end
%     
%     subplot(2,1,1)
%     [F T]=CyclicCorr(ZadoffChu_1,ZadoffChu);
%     h = conj(flipdim(ZadoffChu_1,1));
%     T = conv(ZadoffChu,h);
%     stem(abs(T),'col',colList(i,:))
%     grid on
%     set(gca, 'FontSize',14); 
%     xlabel('Index');
%     ylabel('Magnitude');
%     legend(Leg)
% %     axis([-1 10 0 63])
% 
%     hold on
%     
%     subplot(2,1,2)
%     endpoint = Nrs_zc-shift;
%     temp = ZadoffChu(endpoint+1:end);
%     temp2 = ZadoffChu(1:endpoint);
%     ZC_shift = [temp temp2];
%     [F T]=CyclicCorr(ZadoffChu_1,ZC_shift);
%     stem(0:Nrs_zc-1,abs(T),'col',colList(i,:))
%     grid on
%     set(gca, 'FontSize',14); 
%     xlabel('Index');
%     ylabel('Magnitude');
%     legend(Leg)
% %     axis([-1 10 0 63])
%     hold on
%     
% end 