function [...] = pwelch(x, ...)
  % sort out parameters
  if nargin < 1, 
    usage("[Pxx, w] = pwelch(x,nfft,Fs,window,overlap,pc,range,units,trend)");
  end
  va_start(); 

  % Determine if we are called as pwelch, csd, cohere or tfe
  if isstr(x)
    calledby = x;
  else
    calledby = 'pwelch';
  end
  if !isstr(x)
    ftype = 1;
  elseif strcmp(x, 'csd')
    ftype = 2;
  elseif strcmp(x, 'cohere')
    ftype = 3;
  elseif strcmp(x, 'tfe')
    ftype = 4;
  endif

  ## Sort out x and y vectors
  if ftype!=1, 
    x=va_arg(); y=va_arg(); 
    first = 4;
  else
    y=[];
    first = 2;
  endif
  if (columns(x) != 1 && rows(x) != 1) || ...
    (!isempty(y) && columns(y) != 1 && rows(y) != 1)
    error ([calledby, " data must be a vector"]);
  end
  if columns(x) != 1, x = x'; end
  if columns(y) != 1, y = y'; end
  if !isempty(y) && rows(x)!=rows(y)
    error ([calledby, " x and y vectors must be the same length"]);
  endif

  ## interpret remaining arguments
  trend=nfft=Fs=window=overlap=usewhole=usedB=[];
  ci=0;  ## default no confidence intervals
  pos=0; ## no positional parameters yet interpreted.
  for i=first:nargin
    arg = va_arg();
    if isstr(arg), 
      arg=tolower(arg); 
      if strcmp(arg, 'squared')
      	usedB = 0;
      elseif strcmp(arg, 'db')
	usedB = 1;
      elseif strcmp(arg, 'whole')
	usewhole = 1;
      elseif strcmp(arg, 'half')
	usewhole = 0;
      elseif strcmp(arg, 'none')
      	trend = -1;
      elseif strcmp(arg, 'mean')
      	trend = 0;
      elseif strcmp(arg, 'linear')
      	trend = 1;
      else
      	error([calledby, " doesn't understand '", arg, "'"]);
      endif
    elseif pos == 0
      nfft = arg;
      pos++;
    elseif pos == 1
      Fs = arg;
      pos++;
    elseif pos == 2
      window = arg;
      pos++;
    elseif pos == 3
      overlap = arg;
      pos++;
    elseif pos == 4
      ci = arg;
      pos++;
    else
      usage(usagestr);
    endif
  endfor

  ## Fill in defaults for arguments that aren't specified
  if isempty(nfft), nfft = min(256, length(x)); endif
  if isempty(Fs), Fs = 2; endif
  if isempty(window), window = hanning(nfft); endif
  if isempty(overlap), overlap = length(window)/2; endif
  if isempty(usewhole), usewhole = !isreal(x)||(!isempty(y)&&!isreal(y)); endif
  if isempty(trend), trend=-1; endif
  if isempty(usedB), usedB=ftype!=3; endif # don't default to db for cohere.
  if isempty(ci), ci=0.95; endif # if ci was unspecified, it would be 0

  ## if only the window length is given, generate hanning window
  if length(window) == 1, window = hanning(window); endif
  if rows(window)==1, window = window.'; endif

  ## compute window offsets
  win_size = length(window);
  if (win_size > nfft)
    nfft = win_size;
    warning (sprintf("%s fft size adjusted to %d", calledby, n));
  end
  step = win_size - overlap;

  ## Determine which correlations to compute
  Pci = Pxx = Pyy = Pxy = [];
  if ftype!=2, Pxx = zeros(nfft,1); endif # Not needed for csd
  if ftype==3, Pyy = zeros(nfft,1); endif # Only needed for cohere
  if ftype!=1, Pxy = zeros(nfft,1); endif # Not needed for pwelch
  if ci>0, Pci = zeros(nfft,1); endif     # confidence intervals?

  ## Average the slices
  offset = 1:step:length(x)-win_size+1;
  N = length(offset);
  for i=1:N
    a=x(offset(i):offset(i)+win_size-1);
    if trend>=0, a=detrend(a,trend); endif
    a=fft(postpad(a.*window, nfft));
    if !isempty(Pxx)
      P = a.*conj(a);
      Pxx=Pxx+P;
      if !isempty(Pci), Pci = Pci + P.^2; endif
    endif
    if !isempty(Pxy)
      b=y(offset(i):offset(i)+win_size-1);
      if trend>=0, b=detrend(b,trend); endif
      b=fft(postpad(b.*window, nfft));
      P = a.*conj(b);
      Pxy=Pxy+P;
      if !isempty(Pci), Pci = Pci + P.*conj(P); endif
    endif
    if !isempty(Pyy), Pyy = Pyy + b.*conj(b); endif
  endfor


  if ftype==1,     # pwelch
    P = Pxx;
  elseif ftype==2, # csd
    P = Pxy; 
  elseif ftype==3, # cohere
    P = Pxy.*conj(Pxy)./Pxx./Pyy;
  else             # tfe
    P = Pxy./Pxx;
  endif

  if ftype<=2      # pwelch and csd need normalization and c.i.
    renorm=1/norm(window)^2;
    if ftype==1, renorm=renorm/Fs; endif # only pwelch normalizes frequency
    ## c.i. = mean +/- dev
    ## dev = z_ci*std/sqrt(n)
    ## std = sqrt((N*sum(P^2)-sum(P)^2)/(N*(N-1)))
    ## z_ci = normal_inv( 1-(1-ci)/2 ) = normal_inv( (1+ci)/2 );
    ## normal_inv(x) = sqrt(2) * erfinv(2*x-1)
    ##    => z_ci = sqrt(2)*erfinv(2*(1+ci)/2-1) = sqrt(2)*erfinv 
    ## combining, gives dev as follows:
    if ci>0
      if N>1
      	Pci = sqrt(2*(N*Pci - P.^2)/(N*N*(N-1)))*erfinv(ci);
      else
      	Pci = zeros(nfft,1);
      endif
      Pci=Pci*renorm;
      P = P*renorm/N;
      Pci=[P-Pci, P+Pci];
    else
      P = P*renorm/N;
    endif
  endif
    
  if usedB, 
    P = 10.0*log10(P); 
    Pci = 10.0*log10(Pci);
  endif

  ## extract the positive frequency components
  if usewhole,
    ret_n = nfft;
  elseif rem(nfft,2)==1
    ret_n = (nfft+1)/2;
  else
    ret_n = nfft/2;
  end
  P = P(1:ret_n, :);
  if ci>0,  Pci = Pci(1:ret_n, :); endif
  f = [0:ret_n-1]*Fs/nfft;

  if nargout==0, 
    unwind_protect
      if Fs==2
      	xlabel("Frequency (rad/pi)");
      else
      	xlabel("Frequency (Hz)");
      endif
      if ftype==1
      	title ("Welch's Spectral Estimate Pxx/Fs");
      	ytext="Power Spectral Density";
      elseif ftype==2
      	title ("Cross Spectral Estimate Pxy");
      	ytext="Cross Spectral Density";
      elseif ftype==3
      	title ("Coherence Function Estimate |Pxy|^2/(PxxPyy)");
      	ytext="Coherence ";
      else
      	title ("Transfer Function Estimate Pxy/Pxx");
      	ytext="Transfer";
      endif
      if usedB,
      	ylabel(strcat(ytext, " (dB)"));
      else
      	ylabel(ytext);
      endif
      grid("on");
      if ci>0
      	plot(f, [P, Pci], ";;"); 
      else
      	plot(f, P, ";;");
      endif
    unwind_protect_cleanup
      grid("off");
      title("");
      xlabel("");
      ylabel("");
    end_unwind_protect
  endif
  if nargout>=1, vr_val(P); endif
  if nargout>=2 && ci>0, vr_val(Pci); endif
  if nargout>=2 && ci==0, vr_val(f); endif
  if nargout>=3 && ci>0, vr_val(f); endif

endfunction

