function [ data ] = cellSystem_Config(  )
% Cell Radius
close all; clear all;
model.config.mindataSetSize = 1200;
model.sector_info.R = 500;
model.config.options.BW = [5e6 10e6 15e6 20e6];
model.config.options.BWac = [4.5e6 9e6 13.5e6 18e6];
model.config.options.BSpwr = [43 46 46 46]; %dBm
model.config.options.NumRB = [ 4 8 12 16]; 
model.config.RbFreq = 375e3;
model.config.environment = 'Urban'; % Environment is either urban or rural
model.config.idx = 1;% 2 3 4 
model.config.handover = 3; %dB 
model.config.pwrMargin = 3; %dB 
model.config.sector_size = 3; 
model.config.NumCells = 19;  
model.config.Loaded = 12;
model.config.samplesize = (model.config.sector_size*model.config.NumCells*model.config.Loaded);
model.config.Num_UE = 200000; 
model.config.BW = model.config.options.BW(model.config.idx);
model.config.propagationFreq = 2000e6;
model.config.logStandardDev = 10; %dB
model.config.NumRb = model.config.options.NumRB(model.config.idx);
model.config.NumRbpersub = model.config.options.BWac(model.config.idx)/model.config.RbFreq;
model.config.Num_tier = ceil((model.config.NumCells -1)/6);  
model.config.NumSectors = model.config.sector_size*model.config.NumCells;
model.config.reuse = 1;
model.config.reuseD = sqrt(3*model.config.sector_size)*model.sector_info.R;
model.config.j = 0;  
 
% UE      
model.ue.tx_pwr = 24; % 24dBm Transmit Power 
model.ue.g_rx = 0; % 0 dBi Antenna     
model.ue.NFue = 9; %9dB Noise Figure 
 
% BS
model.bs.tx_pwr = model.config.options.BSpwr(model.config.idx);  
model.bs.NFbs = 5; %dB  
model.bs.MaxPwr_perRB = model.bs.tx_pwr/model.config.NumRb;

[ model ] = getCellinfo( model );  
[ model ] = BSradiation_pattern( model );    
[ model ] = cell_coord( model );  
[ model ] = sector_coord( model );     
[ model ] = createCell( model );      
[ model ] = getAlpaha( model );   
[ data ] = bySector( model ); 

savefile = 'testdata.mat';

try 
    tmp = load(savefile);
    % Downlink
    data.downlink.snir = [tmp.data.downlink.snir; data.downlink.snir];
    data.downlink.thru = [tmp.data.downlink.thru; data.downlink.thru];
    data.downlink.thru_agg = [tmp.data.downlink.thru_agg; data.downlink.thru_agg];
    data.downlink.interference = [tmp.data.downlink.interference; data.downlink.interference];
    
    % Uplink
    data.uplink.snir = [tmp.data.uplink.snir; data.uplink.snir];
    data.uplink.thru = [tmp.data.uplink.thru; data.uplink.thru];
    data.uplink.thru_agg = [tmp.data.uplink.thru_agg; data.uplink.thru_agg];
    data.uplink.interference = [tmp.data.uplink.interference; data.uplink.interference];
    save(savefile, 'data')
catch
    
    save(savefile, 'data')
end
generateCDF( data )


end

