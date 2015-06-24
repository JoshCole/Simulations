function [ r_u_v, q ] = GenReferenceSignal( N_RB_sc, m, u, g, GroupHop, Method )
%   ------------------------------------------------------------------------------------------------
%
%   LTE Reference Signal Simulation Configuration
%
%   function [sim_ctrl,UL_static_cfg,UL_dynamic_cfg,UE_static_cfg,UE_dynamic_cfg,...
%       channel,impairments,testvec] = pusch_config()
%
%   outputs :
%       ~ m                 - 1<= m >= maximum resource blocks in uplink
%       ~ u                 - Group Number 0,1,...,29
%       ~ g                 - Cyclic Shift Index g = 0,1,2,....7
%       ~ GroupHop          - Disable or Enable Group Hop 0, 1 Respectively
%       ~ N_RB_sc            - Number of subcarriers in a resource block 6 or 12
%
%   ------------------------------------------------------------------------------------------------
if strcmp(Method,'PUSCH')
    M_RS_sc = m*N_RB_sc;             % M_RS_sc is the length of the Reference signal and is equal to m.*N_RB_sc
elseif strcmp(Method,'SRS')
    M_RS_sc = m*N_RB_sc/2;           % M_RS_sc is the length of the Reference signal and is equal to m.*N_RB_sc/2 uses a minimum of 4 resource blocks, everyother subcarrier is used hence 6 in one RB
end

% u = 4;%28mod(fgh_ns+fss,30);                     % Group Number
% 0,1,...,29

% Group Hop
% ---------

if GroupHop
    c_i = 0;
    for i = 0:7
        c = GoldSequence(8*ns+i).*2^i;
        c_i = c+c_i;
    end
    fgh_ns = mod(ci,30);
else
    fgh_ns = 0;
end
N_Cell_ID = 0;                              % Cell ID 0, 1, 2,...., 503
% PUCCH
% fss = mod(n_RS_ID,30);
% n_RS_ID = N_Cell_ID;                        % If no value for n_PUCCH_ID
% is configure by higher layers
% Otherwise
% n_RS_ID = n_PUCCH_ID

% PUSCH
% fss = mod(N_Cell_ID+Delta_ss,30);

n_RS_ID = N_Cell_ID;                        % If no value for n_PUSCH_ID is configured by higher layers or if the temporary C-RNTI was used 
                                            % to transmit the most recent uplink-related DCI for the transportblock associated with the corresponding PUSH transmission
% Otherwise
% n_RS_ID = n_PUSCH_ID

fss = mod(n_RS_ID,30);
if strcmp(Method,'PUSCH')
    alpha = 2*pi*g/12;                          % PUSH Demodulation Reference Signal
elseif strcmp(Method,'SRS')
    alpha = 2*pi*g/8;                         % Sounding Reference Signal
end

if M_RS_sc >= 6*N_RB_sc
    v=[0,1];                            % Base Sequence within Group Number v = (0,1) if M_RS_sc >= 6*N_RB_sc else v=0
else
    v=0;
end

n = 0:M_RS_sc-1;
if M_RS_sc >= 3*N_RB_sc

    PrimeValues = Prime_2048;               % All Prime numbers < 2048
    N_RS_ZC = find(PrimeValues < M_RS_sc);  % Select all Prime numbers < M_RS_sc
    N_RS_ZC = PrimeValues(N_RS_ZC(end));    % Select the Largest remaining Prime Number

    q_bar = N_RS_ZC.*(u+1)/31;
    q = floor(q_bar+1/2)+v*(-1)^floor(2*q_bar);
    m = mod(n,N_RS_ZC);                 % Is this the cyclic extention providing a better cubic metric than truncation
    x_qm = zeros(length(q),M_RS_sc);
    for numq = 1:length(q)
        x_qm(numq,:) = exp(-1i.*pi.*q(numq).*m.*(m+1)./N_RS_ZC);
    end
    r_bar_u_v = x_qm; 

else
    r_bar_u_v = exp(1i.*PhiTable(u,n,M_RS_sc).*pi/4);
    q = u;
end
r_u_v = zeros(M_RS_sc,length(q));
for numq = 1:length(q)
    r_u_v(:,numq) =exp(1i.*alpha.*n).*r_bar_u_v(numq,:);
end

end

