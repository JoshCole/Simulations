function [ data ] = collectData( model )
% Collect data from the centeral thee sectos i.e. 1,2,3


% [ ThroughputTable ] = Ttable(  );
% snir_ul = round(model.snir);
% snir_dl = round(model.snir_UL);
% % data.alldata.downlink.thru = zeros(1,1);
% for NumSec = 1:3
%     idxLul = ismember(snir_ul,ThroughputTable(:,1) );
%     idxLdl = ismember(snir_dl,ThroughputTable(:,1) );
%     idx_available_dl = find(model.ue.sectorAssign(idxLdl,3)==NumSec);
%     idx_available_ul = find(model.ue.sectorAssign(idxLul,3)==NumSec);
%     Nidxdl = length(idx_available_dl); Nidxul = length(idx_available_ul);
% %     data.alldata.downlink.thru =
% %     padarray(data.alldata.downlink.thru,Nidxdl-1);
% %     data.alldata.downlink.thru = model.DL_thru(idx_available_dl);
%     tmpidxdl = [];
%     while(length(tmpidxdl)<fullLoad)
%         tmpidxdl = unique([tmpidxdl randi([1 Nidxdl],1,fullLoad)]);
%     end
%     tmpidxul = [];
%     while(length(tmpidxul)<fullLoad)
%         tmpidxul = unique([tmpidxul randi([1 Nidxul],1,fullLoad)]);
%     end
%     tmpidxdl = tmpidxdl(1:12);
%     tmpidxul = tmpidxul(1:12);
%     data.perUE.DL_thru(NumSec,:,model.iteration) = model.DL_thru(idx_available_dl(tmpidxdl),1);
%     data.perUE.snir(NumSec,:,model.iteration) = model.snir(idx_available_dl(tmpidxdl));
%     data.perUE.rsrp(NumSec,:,model.iteration) = model.ue.rx_pwr_dB(idx_available_dl(tmpidxdl));
%     
%     data.perUE.UL_thru(NumSec,:,model.iteration) = model.DL_thru(idx_available_ul(tmpidxul),1);
%     data.perUE.snir_ul(NumSec,:,model.iteration) = model.snir(idx_available_ul(tmpidxul));
%     data.perUE.interferance_ul(NumSec,:,model.iteration) = model.interferance(idx_available_ul(tmpidxul));
% end
end

