function [separated_data]=coherent_ML_receiver(received_signal,H,code_name,rate,num_code,modulator)

%HELP: coherent_ML_receiver
%
%Perform Maximum Likelihood Space Time Decoding. The function can be
%computionnaly expensive if the modulation order is too large.
%
%Input: - received signal
%       - channel matrix: H
%       - code_name: The code name ('Alamouti','OSTBC3',...)
%       - code_rate: The code rate ('1','1/2',...)
%       - num_code: The code number (1,2,...)
%       - modulator object
%Ouput: - separated data
%
%Note: This function requires the Matlab Communication Toolbox
%
%Reference: 
%
%[1]E.G. Larsson,P.Stoica. "Space-time block coding for wireless
%communications", Cambridge Press,2003
%
%Programmed by V. Choqueuse (vincent.choqueuse@gmail.com)

%% extract space time block coding information
Rendement=str2num(rate);
[nb_emitters,code_length]=size(space_time_coding(0,code_name,rate,num_code,1));
nb_symbols_block=code_length*str2num(rate);


%% Construction of the Block Alphabet
alphabet=modulator;

%% construct all the symbol combinations
for indice=nb_symbols_block:-1:1
    tmp=repmat(alphabet,length(alphabet)^(nb_symbols_block-indice),length(alphabet)^(indice-1));
    MIMO_alphabet(indice,:)=tmp(:).';
end    
nb_combination=size(MIMO_alphabet,2);
% construct all the block combination
for indice=1:nb_combination
   MIMO_alphabet_temp=MIMO_alphabet(:,indice);
   C(indice,:,:)=space_time_coding(MIMO_alphabet_temp,code_name,rate,num_code);    
end    

%% ML Decoding
%  this receiver is computationnaly expensive for high order modulation and
%  for a large number of symbol per block 

[nb_receiver,N]	= size(received_signal); 
Nb_bloc=N/code_length;
for indice=1:Nb_bloc
    %extrcat one STBC block
    Xv=received_signal(:,(indice-1)*code_length+1:indice*code_length);
    %Minimize the ML metric
    for num_combination=1:nb_combination
        if code_length==1
            C_temp(:,:)=C(num_combination,:).';  %case of spatial multiplexing
        else    
            C_temp(:,:)=C(num_combination,:,:);
        end
        error(num_combination)=norm(Xv-H*C_temp,'fro');
    end
    [mininum_error,index_min]=min(error);
    %keep the best symbol combination
    separated_data(:,indice)=MIMO_alphabet(:,index_min);
end    
