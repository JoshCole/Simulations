function [ model ] = delData( model )
idxDelRow = model.config.idxDel(:,1);
% idxDelRow_Col = model.config.idxDel;

% Delete UEs
model.ue.data.coord(idxDelRow,:) =[];
[model.config.UEt ,~] = size(model.ue.data.coord);
% Delete sector and Cell Assignment Data
% UE
model.ue.sectorAssign(idxDelRow,:) =[];
model.ue.cellAssign(idxDelRow,:) =[];

% Delete UE data
model.ue.data.angle(idxDelRow,:) =[];
model.ue.data.bs_ue_sep(idxDelRow,:) =[];
model.ue.data.pathloss(idxDelRow,:) =[];
model.ue.data.all_rx_pwr_dB(idxDelRow,:) =[];
model.ue.data.Ppusch(idxDelRow,:) =[];
model.ue.data.rx_pwr_dB(idxDelRow,:) =[];
model.ue.data.sc_ue_sep(idxDelRow,:) =[];

% Delete BS data
model.bs.data.rx_pwr_dB(idxDelRow,:) =[];
model.bs.data.all_rx_pwr_dB(idxDelRow,:) =[];
model.bs.data.g_tx(idxDelRow,:) =[];
end

