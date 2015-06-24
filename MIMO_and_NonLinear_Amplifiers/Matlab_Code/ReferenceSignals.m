% close all
% clear all

m = [29 34 25];
Nrs_zc = 63;
colList = hsv(length(m));
figure
hold on
for i = 1:length(m)
    for n = 0:Nrs_zc-1
        ZadoffChu(n+1) = exp(-1i*(pi*m(i)*n*(n+1))/Nrs_zc);
        Leg{i} = ['Root ', num2str(m(i))];
    end
%     h = conj(flipdim(ZadoffChu,1));
%     y = conv(ZadoffChu,h);
    [F T]=CyclicCorr(ZadoffChu,ZadoffChu);
    stem(real(y),'col',colList(i,:))
    
end 
grid on
set(gca, 'FontSize',14); 
xlabel('Index');
ylabel('Amplitude');
legend(Leg)
% 
% hold off
% % M Sequence
% G = [1 0 0 0 0 1 1]; 
% GF = G(2:end); 
% N = length(G)-1;  
% M = 2^N-1; 
% Reg = zeros(1,N); 
% x = zeros(M,1);
% 
% % Initalise Reg 
% Reg(end)=1;   
%  
% for n = 1:M
%     buffer = mod(sum(GF.*Reg),2);      % Caluculate New Bit
%     Reg = [buffer,Reg(1:N-1)];         % Shift Register   
%     x(n) = Reg(N);
% end
% ReferenceSignal = (x.*2-1).';
% [ ReferenceSignalFreq, ReferenceSignalTime ] = CyclicCorr( ReferenceSignal,ReferenceSignal );
% figure
% stem(real(ReferenceSignalTime))
% grid on
% set(gca, 'FontSize',14); 
% xlabel('Index');
% ylabel('Amplitude');
% legend(['M-sequence, M = ',num2str(M)])
% 
% 
% % Group Hop
% % ---------
% GroupHop = 0;                               % Disable or Enable Group Hop 0, 1 Respectively
% 
% if GroupHop
%     c_i = 0;
%     for i = 0:7
%         c = GoldSequence(8*ns+i).*2^i;
%         c_i = c+c_i;
%     end
%     fgh_ns = mod(ci,30);
% else
%     fgh_ns = 0;
% end
% N_Cell_ID = 0;                              % Cell ID 0, 1, 2,...., 503
% % PUCCH
% % fss = mod(n_RS_ID,30);
% % n_RS_ID = N_Cell_ID;                        % If no value for n_PUCCH_ID
% % is configure by higher layers
% % Otherwise
% % n_RS_ID = n_PUCCH_ID
% 
% % PUSCH
% % fss = mod(N_Cell_ID+Delta_ss,30);
% 
% n_RS_ID = N_Cell_ID;                        % If no value for n_PUSCH_ID is configured by higher layers or if the temporary C-RNTI was used 
%                                             % to transmit the most recent uplink-related DCI for the transportblock associated with the corresponding PUSH transmission
% % Otherwise
% % n_RS_ID = n_PUSCH_ID
% 
% fss = mod(n_RS_ID,30);
% 
% N_RB_sc = 12;                               % Number of subcarriers in a resource block
% m = [1, 2, 3, 4, 5 ];                     % 1<= m >= maximum resource blocks in uplink
% M_RS_sc_Array = m.*N_RB_sc;                 % Length of Reference Signals
% Num_M_RS_sc = length(M_RS_sc_Array);        % Number of Reference Signals
% u = 4;%28mod(fgh_ns+fss,30);                     % Group Number 0,1,...,29
% g = 4;%10                                      % Cyclic Shift Index g = 0,1,2,....7
% 
% 
% 
% alpha = 2*pi*g/12;                          % PUSH Demodulation Reference Signal
% % alpha = 2*pi*g/8;                         % Sounding Reference Signal
% 
% figure
% hold on
% colList = hsv(Num_M_RS_sc);
% for i = 1:Num_M_RS_sc
% 
%     M_RS_sc = M_RS_sc_Array(i);             % M_RS_sc is the length of the Reference signal and is equal to m.*N_RB_sc
%     
%     if M_RS_sc >= 6*N_RB_sc
%         v=[0,1];                            % Base Sequence within Group Number v = (0,1) if M_RS_sc >= 6*N_RB_sc else v=0
%     else
%         v=0;
%     end
% 
%     n = 0:M_RS_sc-1;
%     if M_RS_sc >= 3*N_RB_sc
%        
%         N_RS_ZC = Prime_2048;               % All Prime numbers < 2048
%         N_RS_ZC = find(N_RS_ZC < M_RS_sc);  % Select all Prime numbers < M_RS_sc
%         N_RS_ZC = N_RS_ZC(end);             % Select the Largest remaining Prime Number
%         
%         q_bar = N_RS_ZC.*(u+1/31);
%         q = floor(q_bar+1/2)+v.*(-1)^floor(2*q_bar);
%         m = mod(n,N_RS_ZC);                 % Is this the cyclic extention providing a better cubic metric than truncation
%         x_qm = zeros(length(q),M_RS_sc);
%         for numq = 1:length(q)
%             x_qm(numq,:) = exp(-1i*(pi.*q(numq).*m.*(m+1))/N_RS_ZC);
%         end
%         r_bar_u_v = x_qm; 
%         
%     else
%         r_bar_u_v = exp(1i.*PhiTable(u,n,M_RS_sc).*pi/4);
%         q = u;
%     end
%     r_u_v = zeros(M_RS_sc,length(q));
%     for numq = 1:length(q)
%         r_u_v(:,numq) =exp(1i.*alpha.*n).*r_bar_u_v(numq,:);
%     end
% 
% end
% %     [F T]=CyclicCorr(r_u_v,r_u_v);
%     Leg{i} = ['Msc ', num2str(M_RS_sc)];
%     stem(real(r_u_v),'col',colList(i,:))
% grid on
% set(gca, 'FontSize',14); 
% xlabel('Index');
% ylabel('Amplitude');
% legend(Leg)
% hold off
% 
