function h = gaussfir(BT,NT,OF) 
%   GAUSSFIR   Gaussian FIR Pulse-Shaping Filter Design. 
%   H=GAUSSFIR(BT) designs a low pass FIR gaussian pulse-shaping filter. 
%   BT is the 3-dB bandwidth-symbol time product where B is the two-sided 
%   bandwidth in Hertz and T is in seconds. 
% 
%   H=GAUSSFIR(BT,NT) NT is the number of symbol periods between the start 
%   of the filter impulse response and its peak. If NT is not specified,  
%   NT = 3 is used. 
% 
%   H=GAUSSFIR(BT,NT,OF) OF is the oversampling factor, that is, the number 
%   of samples per symbol. If OF is not specified, OF = 2 is used. 
% 
%   The length of the impulse response of the filter is given by 2*OF*NT+1. 
%   Also, the coefficients H are normalized so that the nominal passband 
%   gain is always equal to one. 
% 
%   % EXAMPLE: Design a Gaussian filter to be used in a GSM GMSK scheme. 
%   BT = .3; % 3-dB bandwidth-symbol time 
%   OF = 8;  % Oversampling factor (i.e., number of samples per symbol) 
%   NT = 2;  % 2 symbol periods to the filters peak.  
%   h = gaussfir(BT,NT,OF);  
%   hfvt = fvtool(h,'impulse'); 
% 
%   See also FIRRCOS. 
 
%   References: 
%   [1] Rappaport T.S., Wireless Communications Principles and Practice,   
%   Prentice Hall, 1996 
%   [2] Krishnapura N., Pavan S., Mathiazhagan C., Ramamurthi B., "A 
%   Baseband Pulse Shaping Filter for Gaussian Minimum Shift Keying," 
%   Proceedings of the 1998 IEEE International Symposium on Circuits and 
%   Systems, 1998. ISCAS '98.  
 
%   Author: P. Costa 
%   Copyright 2004 The MathWorks, Inc. 
%   $Revision: 1.1.6.1 $  $Date: 2004/12/26 22:15:57 $ 
 
% Validate number I/O arguments. 
error(nargchk(1,3,nargin)); 
error(nargoutchk(0,1,nargout)); 
 
if nargin < 2, NT = 3; end 
if nargin < 3, OF = 2; end 
 
msg = ''; 
 
% Check for valid BT 
msg = chkBT(BT,'invalidBT','3-dB bandwidth-symbol time product'); 
error(msg) 
 
% Convert to t in which to compute the filter coefficients 
t= convert2t(OF,NT); 
 
% Equation 5.53 of [1] 
alpha = sqrt(2*log(2))/(BT); 
 
% Equation 5.54 of [1] 
h = (sqrt(pi)/alpha)*exp(-(pi^2/alpha^2)*t.^2);  
  
% Normalize coefficients 
h = h./sum(h); 
 
 
%-------------------------------------------------------------------------- 
function t = convert2t(OF,NT) 
 
% Check for valid OF and NT 
msg = chkInput(OF,'invalidOSFactor','Oversampling factor'); 
error(msg) 
 
msg = chkInput(NT,'invalidNumSPeriods','Number of symbol periods'); 
error(msg) 
 
% Filter Length 
filtLen = 2*OF*NT+1; 
t = linspace(-NT,NT,filtLen); 
 
 
%----------------------------------------------------------------------- 
function msg = chkBT(val,id,param) 
 
msg = ''; 
if isempty(val) || length(val) > 1 || ~isa(val,'double') || ... 
        ~isreal(val) || val<=0, 
    msg.identifier = generatemsgid(id); 
    msg.message = [param,' must be a real, positive scalar.']; 
    return; 
end 
 
%----------------------------------------------------------------------- 
function msg = chkInput(val,id,param) 
 
msg = ''; 
if isempty(val) || length(val) > 1 || ~isa(val,'double') || ... 
        ~isreal(val) || val~=round(val) || val<=0, 
    msg.identifier = generatemsgid(id); 
    msg.message = [param,' must be a real, positive integer.']; 
    return; 
end 
 
 
% [EOF] 