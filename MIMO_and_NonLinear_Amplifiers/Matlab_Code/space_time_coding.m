function [modulated_coded_symbols]=space_time_coding(modulated_symbols,code_name,rate,num_code,zero_padding)

%% Space_time_coding
%
%Encode an input stream into a specified space time block code.
%
%Input:     - modulated_symbols: a vector of symbols (if the length does not match with the STBC
%           structure you must use the zero_padding option)
%           - code_name: A string corresponding to the code name
%           ('SM','Alamouti','OSTBC3','STBC3','QOSTBC4','OSTBC4','STBC4','OSTBC5','OSTBC6','OSTBC7')
%           - rate: a string with the code rate '1/2','3/4','1','N'
%           - zero_padding (optionnal): if there is a 5th input and if the length does not
%           match with the space time bloc coding structure, than the modulated symbols are padded with zero. 
%
%Output:    modulated_coded_symbols
%           
%Example: 
%>> C=space_time_coding([1+i,2+2*i],'Alamouti','1',1)
%C =
%   1.0000 + 1.0000i  -2.0000 + 2.0000i
%   2.0000 + 2.0000i   1.0000 - 1.0000i
%
%Remarks: 
%To obtain the code matrix, enter for example:
%>>space_time_coding([1:4]+i*[1:4],'OSTBC4','3/4',1)
%
%To obtain the code length and the number of transmitter used by a STBC,
%enter: size(space_time_coding(0,'OSTBC4','3/4',1,1))
%
%Reference: 
%
%[1] S.M Alamouti " A simple transmitter diversity scheme for wireless
%communications" IEEE J.Select Areas Communication vol 16,
%October 1998
%[2] V.Tarokh "Space Time BLock Codes from Orthogonal Designs" IEEE
%Transactions on Information theory, vol 45, July 1999
%[3] G.Ganesan and P. Stoica "Space time Block Codes: A Maximum SNR
%approach", IEEE Transactions on Information Theory, vol 47, may 2001.
%[4] I. Jafarkhani: "A Quasi Orthogonal Space Time Block Code"
%IEEE Transactions Letters on Communications, Vol. 49,
%No.1, 2001
%
%
%Programmed by V. Choqueuse (vincent.choqueuse@gmail.com)

modulated_coded_symbols=[];

switch code_name

    case 'SM' 
        %% Spatial Multiplexing %%
        nb_emitters=str2num(rate);
        
        if(mod(length(modulated_symbols),nb_emitters)~=0) && (nargin==5)
           modulated_symbols(end+1:end+(nb_emitters-mod(length(modulated_symbols),nb_emitters)))=0;
        end   
        
        if (mod(nb_emitters,1)==0)
            modulated_coded_symbols=reshape(modulated_symbols,nb_emitters,prod(size(modulated_symbols))/nb_emitters);
        else
            fprintf('Error: Spatial Mulitplexing -> code rate must be an integer\n');
        end   
        
    case 'Alamouti'  

        if (strcmp(rate,'1')==1) 
           
            %% Alamouti code %%
            
            %[1] S.M Alamouti " A simple transmitter diversity scheme for wireless
            %communications" IEEE J.Select Areas Communication vol 16, October 1998
            if(mod(length(modulated_symbols),2)~=0) && (nargin==5)
                modulated_symbols(end+1:end+(2-mod(length(modulated_symbols),2)))=0;
            end  
         
            symboles_reshape=reshape(modulated_symbols,2,length(modulated_symbols)/2);
            x1=[1 0];
            x2=[0 1];
            P1=[x1;x2];
            P2=[-x2;x1];
            P= [P1;conj(P2)];
            modulated_coded_symbols=P*symboles_reshape;
            modulated_coded_symbols(3:4,:)=conj(modulated_coded_symbols(3:4,:));
            modulated_coded_symbols=reshape(modulated_coded_symbols,2,prod(size(modulated_coded_symbols))/2);
            
        else
            fprintf('Error: Alamouti coding -> code rate must be 1\n');
        end    
    
    case 'OSTBC3'  
        switch rate
            case '1/2' 
                %% Tarokh OSTBC
                
                %[2] V.Tarokh "Space Time BLock Codes from Orthogonal Designs" IEEE
                %Transactions on Information theory, vol 45, July 1999
                if(mod(length(modulated_symbols),4)~=0) && (nargin==5)
                    modulated_symbols(end+1:end+(4-mod(length(modulated_symbols),4)))=0;
                end   
                
                symboles_reshape=reshape(modulated_symbols,4,length(modulated_symbols)/4);
                x1=[1 0 0 0];
                x2=[0 1 0 0];
                x3=[0 0 1 0];
                x4=[0 0 0 1];
                P1=[x1;x2; x3];
                P2=[-x2;x1;-x4];
                P3=[-x3;x4;x1];
                P4=[-x4;-x3;x2];
                P=[P1;P2;P3;P4];
                modulated_coded_symbols=[P;P]*symboles_reshape;
                modulated_coded_symbols(13:end,:)=conj(modulated_coded_symbols(13:end,:));
                modulated_coded_symbols=reshape(modulated_coded_symbols,3,prod(size(modulated_coded_symbols))/3);
                
            case '3/4' 
                switch num_code
                    case 1
                        
                        %% Ganesan STBC
                        
                        %[3] G.Ganesan and P. Stoica "Space time Block Codes: A Maximum SNR
                        %approach", IEEE Transactions on Information Theory, vol 47, may 2001.
                        
                        if(mod(length(modulated_symbols),3)~=0) && (nargin==5)
                            modulated_symbols(end+1:end+(3-mod(length(modulated_symbols),3)))=0;
                        end   
                        
                        symboles_reshape=reshape(modulated_symbols,3,length(modulated_symbols)/3);
                        symboles_reshape_conj=conj(symboles_reshape);
                        x1=[1 0 0];
                        x2=[0 1 0];
                        x3=[0 0 1];
                        x0=[0 0 0];
                
                        P1=[x1;x0;x0];
                        P2=[x0;x1;-x3];
                        P3=[x2;x0;x0];
                        P4=[-x3;x0;x0];
                        P=[P1;P2;P3;P4];
                
                        P1_conj=[x0;x0;-x2];
                        P2_conj=[x0;x0;x0];
                        P3_conj=[x0;x3;x1];
                        P4_conj=[x0;x2;x0];
                        P_conj=[P1_conj;P2_conj;P3_conj;P4_conj];
                
                        modulated_coded_symbols=P*symboles_reshape+P_conj*symboles_reshape_conj;
                        modulated_coded_symbols=reshape(modulated_coded_symbols,3,prod(size(modulated_coded_symbols))/3);
                    case 2
                        
                        %% Ganesan STBC
                        
                        %[3] G.Ganesan and P. Stoica "Space time Block Codes: A Maximum SNR
                        %approach", IEEE Transactions on Information
                        %Theory, vol 47, may 2001.
                        
                        if(mod(length(modulated_symbols),3)~=0) && (nargin==5)
                            modulated_symbols(end+1:end+(3-mod(length(modulated_symbols),3)))=0;
                        end  
                        
                        symboles_reshape=reshape(modulated_symbols,3,length(modulated_symbols)/3);
                        symboles_reshape_conj=conj(symboles_reshape);
                        x1=[1 0 0];
                        x2=[0 1 0];
                        x3=[0 0 1];
                        x0=[0 0 0];
                
                        P1=[x1;x2;x3];
                        P2=[x0;x0;x0];
                        P3=[x0;x0;x0];
                        P4=[x0;x0;x0];
                        P=[P1;P2;P3;P4];
                
                        P1_conj=[x0;x0;x0];
                        P2_conj=[-x2;x1;x0];
                        P3_conj=[x3;x0;-x1];
                        P4_conj=[x0;-x3;x2];
                        P_conj=[P1_conj;P2_conj;P3_conj;P4_conj];
                
                        modulated_coded_symbols=P*symboles_reshape+P_conj*symboles_reshape_conj;
                        modulated_coded_symbols=reshape(modulated_coded_symbols,3,prod(size(modulated_coded_symbols))/3);
                    case 3   
                         %% Tarokh OSTBC
                
                        %[2] V.Tarokh "Space Time BLock Codes from Orthogonal Designs" IEEE
                        %Transactions on Information theory, vol 45, July 1999
                        
                        if(mod(length(modulated_symbols),3)~=0) && (nargin==5)
                            modulated_symbols(end+1:end+(3-mod(length(modulated_symbols),3)))=0;
                        end  
                        
                        symboles_reshape=reshape(modulated_symbols,3,length(modulated_symbols)/3);
                        symboles_reshape_conj=conj(symboles_reshape);
                        x1=[1 0 0];
                        x2=[0 1 0];
                        x3=[0 0 1];
                        x0=[0 0 0];
                
                        P1=[x1;x2;x3/sqrt(2)];
                        P2=[x0;x0;x3/sqrt(2)];
                        P3=[x0;x0;(-x1+x2)/2];
                        P4=[x0;x0;(x2+x1)/2];
                        P=[P1;P2;P3;P4];
                
                        P1_conj=[x0;x0;x0];
                        P2_conj=[-x2;x1;x0];
                        P3_conj=[x3/sqrt(2);x3/sqrt(2);(-x1-x2)/2];
                        P4_conj=[x3/sqrt(2);-x3/sqrt(2);(x2-x1)/2];
                        P_conj=[P1_conj;P2_conj;P3_conj;P4_conj];
                
                        modulated_coded_symbols=P*symboles_reshape+P_conj*symboles_reshape_conj;
                        modulated_coded_symbols=reshape(modulated_coded_symbols,3,prod(size(modulated_coded_symbols))/3);
                    otherwise
                       fprintf(sprintf('Error: the %i OSTBC3 with rate %s is not currently implemented\n',num_code,rate)); 
                     
                end    
                          
            otherwise
                fprintf(sprintf('Error: OSTBC3 with rate %s is not currently implemented\n',rate));
        end       
     case 'STBC3' 
        switch rate 
            case '2'   
                
               %% IEEE 802.16e (non orthogonal) STBC 
                
               if(mod(length(modulated_symbols),4)~=0) && (nargin==5)
                     modulated_symbols(end+1:end+(4-mod(length(modulated_symbols),4)))=0;
               end  
                
               symboles_reshape=reshape(modulated_symbols,4,length(modulated_symbols)/4);
               symboles_reshape_conj=conj(symboles_reshape);
               x1=[1 0 0 0];
               x2=[0 1 0 0];
               x3=[0 0 1 0];
               x4=[0 0 0 1];
               x0=[0 0 0 0];
               
               P1=[x1;x0;x3];
               P2=[x2;x0;x4];
               P=[P1;P2];
               P1_conj=[x0;-x2;x0];
               P2_conj=[x0;x1;x0];
               P_conj=[P1_conj;P2_conj];
                        
               modulated_coded_symbols=P*symboles_reshape+P_conj*symboles_reshape_conj;
               modulated_coded_symbols=reshape(modulated_coded_symbols,3,prod(size(modulated_coded_symbols))/3);
           otherwise
                fprintf(sprintf('Error: STBC3 with rate %s is not currently implemented\n',rate));      
        end
                          
    case 'QOSTBC4'
       
        switch rate 
            case '1'  
                %% Quasi Orthogonal STBC
        
                %[4] I. Jafarkhani: "A Quasi Orthogonal Space Time Block Code"
                %IEEE Transactions Letters on Communications, Vol. 49, No.1, 2001
                switch num_code
                    case 1
                        
                        if(mod(length(modulated_symbols),4)~=0) && (nargin==5)
                            modulated_symbols(end+1:end+(4-mod(length(modulated_symbols),4)))=0;
                        end  
                        
                        symboles_reshape=reshape(modulated_symbols,4,length(modulated_symbols)/4);
                        symboles_reshape_conj=conj(symboles_reshape);
                        x1=[1 0 0 0];
                        x2=[0 1 0 0];
                        x3=[0 0 1 0];
                        x4=[0 0 0 1];
                        x0=[0 0 0 0];
                
                        P1=[x1;x0;x0;x4];
                        P2=[x2;x0;x0;-x3];
                        P3=[x3;x0;x0;-x2];
                        P4=[x4;x0;x0;x1];
                        P=[P1;P2;P3;P4];
                
                        P1_conj=[x0;-x2;-x3;x0];
                        P2_conj=[x0;x1;-x4;x0];
                        P3_conj=[x0;-x4;x1;x0];
                        P4_conj=[x0;x3;x2;x0];
                        P_conj=[P1_conj;P2_conj;P3_conj;P4_conj];
                
                        modulated_coded_symbols=P*symboles_reshape+P_conj*symboles_reshape_conj;
                        modulated_coded_symbols=reshape(modulated_coded_symbols,4,prod(size(modulated_coded_symbols))/4);
                end    
            otherwise
                fprintf(sprintf('Error: QOSTBC4 with rate %s does not exits',rate));
        end       
    
    case 'OSTBC4'
        switch rate 
            case '1/2'     
                %% Tarokh OSTBC
                
                %[2] V.Tarokh "Space Time BLock Codes from Orthogonal Designs" IEEE
                %Transactions on Information theory, vol 45, July 1999
                
                if(mod(length(modulated_symbols),4)~=0) && (nargin==5)
                     modulated_symbols(end+1:end+(4-mod(length(modulated_symbols),4)))=0;
                end  
                
                symboles_reshape=reshape(modulated_symbols,4,length(modulated_symbols)/4);
                x1=[1 0 0 0];
                x2=[0 1 0 0];
                x3=[0 0 1 0];
                x4=[0 0 0 1];
                
                P1=[x1;x2;x3;x4];
                P2=[-x2;x1;-x4;x3];
                P3=[-x3;x4;x1;-x2];
                P4=[-x4;-x3;x2;x1];
                P=[P1;P2;P3;P4];
                modulated_coded_symbols=[P;P]*symboles_reshape;
                modulated_coded_symbols(17:end,:)=conj(modulated_coded_symbols(17:end,:));
                modulated_coded_symbols=reshape(modulated_coded_symbols,4,prod(size(modulated_coded_symbols))/4);
            
            case '3/4'  
                switch num_code
                    case 1
                        
                        %% Ganesan STBC
                        
                        %[3] G.Ganesan and P. Stoica "Space time Block Codes: A Maximum SNR
                        %approach", IEEE Transactions on Information
                        %Theory, vol 47, may 2001.
                        
                        if(mod(length(modulated_symbols),3)~=0) && (nargin==5)
                            modulated_symbols(end+1:end+(3-mod(length(modulated_symbols),3)))=0;
                        end  
                        
                        symboles_reshape=reshape(modulated_symbols,3,length(modulated_symbols)/3);
                        symboles_reshape_conjuge=conj(symboles_reshape);
                        x1=[1 0 0];
                        x2=[0 1 0];
                        x3=[0 0 1];
                        x0=[0 0 0];

                        P1=[x1;x0;x0;x0];
                        P2=[x0;x1;-x3;-x2];
                        P3=[x2;x0;x0;x0];
                        P4=[-x3;x0;x0;x0];
                        P=[P1;P2;P3;P4];
                
                        P1_conj=[x0;x0;-x2;x3];
                        P2_conj=[x0;x0;x0;x0];
                        P3_conj=[x0;x3;x1;x0];
                        P4_conj=[x0;x2;x0;x1];
                        P_conj=[P1_conj;P2_conj;P3_conj;P4_conj];
                        modulated_coded_symbols=P*symboles_reshape+P_conj*symboles_reshape_conjuge;
                        modulated_coded_symbols=reshape(modulated_coded_symbols,4,prod(size(modulated_coded_symbols))/4);
                        
                    case 2 
                        
                         %% Tarokh OSTBC
                
                         %[2] V.Tarokh "Space Time BLock Codes from Orthogonal Designs" IEEE
                         %Transactions on Information theory, vol 45, July 1999
                        
                        if(mod(length(modulated_symbols),3)~=0) && (nargin==5)
                            modulated_symbols(end+1:end+(3-mod(length(modulated_symbols),3)))=0;
                        end  
                        
                        symboles_reshape=reshape(modulated_symbols,3,length(modulated_symbols)/3);
                        symboles_reshape_conj=conj(symboles_reshape);
                        x1=[1 0 0];
                        x2=[0 1 0];
                        x3=[0 0 1];
                        x0=[0 0 0];
                
                        P1=[x1;x2;x3/sqrt(2);x3/sqrt(2)];
                        P2=[x0;x0;x3/sqrt(2);-x3/sqrt(2)];
                        P3=[x0;x0;(-x1+x2)/2;(-x2+x1)/2];
                        P4=[x0;x0;(x2+x1)/2;-(x1+x2)/2];
                        P=[P1;P2;P3;P4];
                
                        P1_conj=[x0;x0;x0;x0];
                        P2_conj=[-x2;x1;x0;x0];
                        P3_conj=[x3/sqrt(2);x3/sqrt(2);(-x1-x2)/2;(-x2-x1)/2];
                        P4_conj=[x3/sqrt(2);-x3/sqrt(2);(x2-x1)/2;(-x1+x2)/2];
                        P_conj=[P1_conj;P2_conj;P3_conj;P4_conj];
                
                        modulated_coded_symbols=P*symboles_reshape+P_conj*symboles_reshape_conj;
                        modulated_coded_symbols=reshape(modulated_coded_symbols,4,prod(size(modulated_coded_symbols))/4);
                        
                    otherwise
                       fprintf(sprintf('Error: the %i OSTBC3 with rate %s is not currently implemented\n',num_code,rate)); 
                        
                end
            otherwise
                fprintf(sprintf('Error: OSTBC4 with rate %s is not currently implemented\n',rate));      
        end
    %----------------------%
    %         STBC4        %
    %----------------------% 
    case 'STBC4' 
        switch rate 
           case '1'    
                switch num_code
                    case 1
                        
                        %% IEEE 802.16e (non orthogonal) STBC 
                        
                        if(mod(length(modulated_symbols),4)~=0) && (nargin==5)
                            modulated_symbols(end+1:end+(4-mod(length(modulated_symbols),4)))=0;
                        end  
                        
                        symboles_reshape=reshape(modulated_symbols,4,length(modulated_symbols)/4);
                        symboles_reshape_conj=conj(symboles_reshape);
                        x1=[1 0 0 0];
                        x2=[0 1 0 0];
                        x3=[0 0 1 0];
                        x4=[0 0 0 1];
                        x0=[0 0 0 0];
                        
                        P1=[x1;x2;x0;x0];
                        P2=[x0;x0;x0;x0];
                        P3=[x0;x0;x3;x4];
                        P4=[x0;x0;x0;x0];
                        P=[P1;P2;P3;P4];
                        P1_conj=[x0;x0;x0;x0];
                        P2_conj=[-x2;x1;x0;x0];
                        P3_conj=[x0;x0;x0;x0];
                        P4_conj=[x0;x0;-x4;x3];
                        P_conj=[P1_conj;P2_conj;P3_conj;P4_conj];
                        
                        modulated_coded_symbols=P*symboles_reshape+P_conj*symboles_reshape_conj;
                        modulated_coded_symbols=reshape(modulated_coded_symbols,4,prod(size(modulated_coded_symbols))/4);
                        
                    case 2
                        
                        %% IEEE 802.16e (non orthogonal) STBC 
                        
                        if(mod(length(modulated_symbols),4)~=0) && (nargin==5)
                            modulated_symbols(end+1:end+(4-mod(length(modulated_symbols),4)))=0;
                        end  
                        
                        symboles_reshape=reshape(modulated_symbols,4,length(modulated_symbols)/4);
                        symboles_reshape_conj=conj(symboles_reshape);
                        x1=[1 0 0 0];
                        x2=[0 1 0 0];
                        x3=[0 0 1 0];
                        x4=[0 0 0 1];
                        x0=[0 0 0 0];
                        
                        P1=[x1;x0;x2;x0];
                        P2=[x0;x0;x0;x0];
                        P3=[x0;x3;x0;x4];
                        P4=[x0;x0;x0;x0];
                        P=[P1;P2;P3;P4];
                        P1_conj=[x0;x0;x0;x0];
                        P2_conj=[-x2;x0;x1;x0];
                        P3_conj=[x0;x0;x0;x0];
                        P4_conj=[x0;-x4;x0;x3];
                        P_conj=[P1_conj;P2_conj;P3_conj;P4_conj];
                        
                        modulated_coded_symbols=P*symboles_reshape+P_conj*symboles_reshape_conj;
                        modulated_coded_symbols=reshape(modulated_coded_symbols,4,prod(size(modulated_coded_symbols))/4);
                    case 3
                        
                        %% IEEE 802.16e (non orthogonal) STBC 
                        
                        if(mod(length(modulated_symbols),4)~=0) && (nargin==5)
                            modulated_symbols(end+1:end+(4-mod(length(modulated_symbols),4)))=0;
                        end  
                        
                        symboles_reshape=reshape(modulated_symbols,4,length(modulated_symbols)/4);
                        symboles_reshape_conj=conj(symboles_reshape);
                        x1=[1 0 0 0];
                        x2=[0 1 0 0];
                        x3=[0 0 1 0];
                        x4=[0 0 0 1];
                        x0=[0 0 0 0];
                        
                        P1=[x1;x0;x0;x2];
                        P2=[x0;x0;x0;x0];
                        P3=[x0;x3;x4;x0];
                        P4=[x0;x0;x0;x0];
                        P=[P1;P2;P3;P4];
                        P1_conj=[x0;x0;x0;x0];
                        P2_conj=[-x2;x0;x0;x1];
                        P3_conj=[x0;x0;x0;x0];
                        P4_conj=[x0;-x4;-x3;x0];
                        P_conj=[P1_conj;P2_conj;P3_conj;P4_conj];
                        
                        modulated_coded_symbols=P*symboles_reshape+P_conj*symboles_reshape_conj;
                        modulated_coded_symbols=reshape(modulated_coded_symbols,4,prod(size(modulated_coded_symbols))/4);
                    otherwise
                       fprintf(sprintf('Error: the %i STBC4 with rate %s is not currently implemented\n',num_code,rate)); 
                end

             case '2'
                 switch num_code
                    case 1
                        
                        %% IEEE 802.16e (non orthogonal) STBC 
                        
                        if(mod(length(modulated_symbols),8)~=0) && (nargin==5)
                            modulated_symbols(end+1:end+(8-mod(length(modulated_symbols),8)))=0;
                        end  
                        
                        symboles_reshape=reshape(modulated_symbols,8,length(modulated_symbols)/8);
                        symboles_reshape_conj=conj(symboles_reshape);
                        x1=[1 0 0 0 0 0 0 0];
                        x2=[0 1 0 0 0 0 0 0];
                        x3=[0 0 1 0 0 0 0 0];
                        x4=[0 0 0 1 0 0 0 0];
                        x5=[0 0 0 0 1 0 0 0];
                        x6=[0 0 0 0 0 1 0 0];
                        x7=[0 0 0 0 0 0 1 0];
                        x8=[0 0 0 0 0 0 0 1];
                        x0=[0 0 0 0 0 0 0 0];
                        
                        P1=[x1;x2;x3;x4];
                        P2=[x0;x0;x0;x0];
                        P3=[x5;x6;x7;x8];
                        P4=[x0;x0;x0;x0];
                        P=[P1;P2;P3;P4];
                        P1_conj=[x0;x0;x0;x0];
                        P2_conj=[-x2;x1;-x4;x3];
                        P3_conj=[x0;x0;x0;x0];
                        P4_conj=[-x7;-x8;x5;x6];
                        P_conj=[P1_conj;P2_conj;P3_conj;P4_conj];
                        
                        modulated_coded_symbols=P*symboles_reshape+P_conj*symboles_reshape_conj;
                        modulated_coded_symbols=reshape(modulated_coded_symbols,4,prod(size(modulated_coded_symbols))/4);
                  case 2
                        
                        %% IEEE 802.16e (non orthogonal) STBC 
                      
                        if(mod(length(modulated_symbols),4)~=0) && (nargin==5)
                            modulated_symbols(end+1:end+(4-mod(length(modulated_symbols),4)))=0;
                        end  
                      
                        symboles_reshape=reshape(modulated_symbols,4,length(modulated_symbols)/4);
                        symboles_reshape_conj=conj(symboles_reshape);
                        x1=[1 0 0 0 ];
                        x2=[0 1 0 0 ];
                        x3=[0 0 1 0 ];
                        x4=[0 0 0 1 ];
                        x0=[0 0 0 0 ];
                        
                        P1=[x1;x0;x3;x0];
                        P2=[x2;x0;x4;x0];
                        P=[P1;P2];
                        P1_conj=[x0;-x2;x0;-x4];
                        P2_conj=[x0;x1;x0;x3];
                        P_conj=[P1_conj;P2_conj];
                        
                        modulated_coded_symbols=P*symboles_reshape+P_conj*symboles_reshape_conj;
                        modulated_coded_symbols=reshape(modulated_coded_symbols,4,prod(size(modulated_coded_symbols))/4);   
                   otherwise
                       fprintf(sprintf('Error: the %i STBC4 with rate %s is not currently implemented\n',num_code,rate)); 
                 end 
            case '3'
                 switch num_code
                    case 1
                        
                        %% IEEE 802.16e (non orthogonal) STBC 
                        
                        if(mod(length(modulated_symbols),6)~=0) && (nargin==5)
                            modulated_symbols(end+1:end+(6-mod(length(modulated_symbols),6)))=0;
                        end  
                        
                        symboles_reshape=reshape(modulated_symbols,6,length(modulated_symbols)/6);
                        symboles_reshape_conj=conj(symboles_reshape);
                        x1=[1 0 0 0 0 0];
                        x2=[0 1 0 0 0 0];
                        x3=[0 0 1 0 0 0];
                        x4=[0 0 0 1 0 0];
                        x5=[0 0 0 0 1 0];
                        x6=[0 0 0 0 0 1];
                        x0=[0 0 0 0 0 0];
                        
                        P1=[x1;x0;x3;x5];
                        P2=[x2;x0;x4;x6];
                        P=[P1;P2];
                        P1_conj=[x0;-x2;x0;x0];
                        P2_conj=[x0;x1;x0;x0];
                        P_conj=[P1_conj;P2_conj];
                        
                        modulated_coded_symbols=P*symboles_reshape+P_conj*symboles_reshape_conj;
                        modulated_coded_symbols=reshape(modulated_coded_symbols,4,prod(size(modulated_coded_symbols))/4);   
                   otherwise
                       fprintf(sprintf('Error: the %i STBC4 with rate %s is not currently implemented\n',num_code,rate)); 
                 end     
           otherwise
                fprintf(sprintf('Error: STBC4 with rate %s is not currently implemented\n',rate));      
    end
  
    case 'OSTBC5' 
         switch rate 
            case '1/2' 
                 %% Tarokh OSTBC
                
                %[2] V.Tarokh "Space Time BLock Codes from Orthogonal Designs" IEEE
                %Transactions on Information theory, vol 45, July 1999 
                
               if(mod(length(modulated_symbols),8)~=0) && (nargin==5)
                      modulated_symbols(end+1:end+(8-mod(length(modulated_symbols),8)))=0;
               end  
                
                
                symboles_reshape=reshape(modulated_symbols,8,length(modulated_symbols)/8);
                x1=[1 0 0 0 0 0 0 0];
                x2=[0 1 0 0 0 0 0 0];
                x3=[0 0 1 0 0 0 0 0];
                x4=[0 0 0 1 0 0 0 0];
                x5=[0 0 0 0 1 0 0 0];
                x6=[0 0 0 0 0 1 0 0];
                x7=[0 0 0 0 0 0 1 0];
                x8=[0 0 0 0 0 0 0 1];
                P1=[x1;x2;x3;x4;x5];
                P2=[-x2;x1;x4;-x3;x6];
                P3=[-x3;-x4;x1;x2;x7];
                P4=[-x4;x3;-x2;x1;x8];
                P5=[-x5;-x6;-x7;-x8;x1];
                P6=[-x6;x5;-x8;x7;-x2];
                P7=[-x7;x8;x5;-x6;-x3];
                P8=[-x8;-x7;x6;x5;-x4];
                P=[P1;P2;P3;P4;P5;P6;P7;P8];
                modulated_coded_symbols=[P;P]*symboles_reshape;
                modulated_coded_symbols(41:end,:)=conj(modulated_coded_symbols(41:end,:));
                modulated_coded_symbols=reshape(modulated_coded_symbols,5,prod(size(modulated_coded_symbols))/5);
            otherwise
                fprintf(sprintf('Error: OSTBC5 with rate %s is not currently implemented\n',rate));
         end
         
        
    %----------------------%
    %         OSTBC6       %
    %----------------------%      
    case 'OSTBC6' 
         switch rate 
            case '1/2'  
                %% Tarokh OSTBC
                
                %[2] V.Tarokh "Space Time BLock Codes from Orthogonal Designs" IEEE
                %Transactions on Information theory, vol 45, July 1999  
                
                if(mod(length(modulated_symbols),8)~=0) && (nargin==5)
                      modulated_symbols(end+1:end+(8-mod(length(modulated_symbols),8)))=0;
               end  

                symboles_reshape=reshape(modulated_symbols,8,length(modulated_symbols)/8);
                x1=[1 0 0 0 0 0 0 0];
                x2=[0 1 0 0 0 0 0 0];
                x3=[0 0 1 0 0 0 0 0];
                x4=[0 0 0 1 0 0 0 0];
                x5=[0 0 0 0 1 0 0 0];
                x6=[0 0 0 0 0 1 0 0];
                x7=[0 0 0 0 0 0 1 0];
                x8=[0 0 0 0 0 0 0 1];
                P1=[x1;x2;x3;x4;x5;x6];
                P2=[-x2;x1;x4;-x3;x6;-x5];
                P3=[-x3;-x4;x1;x2;x7;x8];
                P4=[-x4;x3;-x2;x1;x8;-x7];
                P5=[-x5;-x6;-x7;-x8;x1;x2];
                P6=[-x6;x5;-x8;x7;-x2;x1];
                P7=[-x7;x8;x5;-x6;-x3;x4];
                P8=[-x8;-x7;x6;x5;-x4;-x3];
                P=[P1;P2;P3;P4;P5;P6;P7;P8];
                modulated_coded_symbols=[P;P]*symboles_reshape;
                modulated_coded_symbols(49:end,:)=conj(modulated_coded_symbols(49:end,:));
                modulated_coded_symbols=reshape(modulated_coded_symbols,6,prod(size(modulated_coded_symbols))/6);
            otherwise
                fprintf(sprintf('Error: OSTBC6 with rate %s is not currently implemented\n',rate));
         end
         
    %----------------------%
    %         OSTBC7       %
    %----------------------%         
    case 'OSTBC7'
        switch rate 
            case '1/2'   
                
                %% Tarokh OSTBC
                
                %[2] V.Tarokh "Space Time BLock Codes from Orthogonal Designs" IEEE
                %Transactions on Information theory, vol 45, July 1999  
                
                if(mod(length(modulated_symbols),8)~=0) && (nargin==5)
                      modulated_symbols(end+1:end+(8-mod(length(modulated_symbols),8)))=0;
                end  
                
                symboles_reshape=reshape(modulated_symbols,8,length(modulated_symbols)/8);
                x1=[1 0 0 0 0 0 0 0];
                x2=[0 1 0 0 0 0 0 0];
                x3=[0 0 1 0 0 0 0 0];
                x4=[0 0 0 1 0 0 0 0];
                x5=[0 0 0 0 1 0 0 0];
                x6=[0 0 0 0 0 1 0 0];
                x7=[0 0 0 0 0 0 1 0];
                x8=[0 0 0 0 0 0 0 1];
                P1=[x1;x2;x3;x4;x5;x6;x7];
                P2=[-x2;x1;x4;-x3;x6;-x5;-x8];
                P3=[-x3;-x4;x1;x2;x7;x8;-x5];
                P4=[-x4;x3;-x2;x1;x8;-x7;x6];
                P5=[-x5;-x6;-x7;-x8;x1;x2;x3];
                P6=[-x6;x5;-x8;x7;-x2;x1;-x4];
                P7=[-x7;x8;x5;-x6;-x3;x4;x1];
                P8=[-x8;-x7;x6;x5;-x4;-x3;x2];
                P=[P1;P2;P3;P4;P5;P6;P7;P8];
                modulated_coded_symbols=[P;P]*symboles_reshape;
                modulated_coded_symbols(57:end,:)=conj(modulated_coded_symbols(57:end,:));
                modulated_coded_symbols=reshape(modulated_coded_symbols,7,prod(size(modulated_coded_symbols))/7);    
            otherwise
                fprintf(sprintf('Error: OSTBC7 with rate %s is not currently implemented\n',rate));
        end
    
    otherwise
        fprintf('Error: Unknown Linear Space Time Bloc Code\n');
        
end       