function in = input_sspp(time,intype,inparam)

totsamp = length(time);

% Spike inputs -----------------------------------------------------------%
if strmatch('spike',intype)
    rows  = inparam(1);
    ratio = inparam(2:end);
    in = zeros(rows,totsamp);
    for row = 1:rows
        in(row,ratio(row):ratio(row):end) = 1;
    end
end
%-------------------------------------------------------------------------%

% Harmonic excitation signals --------------------------------------------%
if strmatch('harmonic',intype)
    num  = length(inparam); 
    in   = zeros(2*num,totsamp);
    for ii = 1:num
        blk = (ii-1)*2+1:ii*2;
        in(blk,:) = [cos(inparam(ii)*time);sin(inparam(ii)*time)];
    end
end
%-------------------------------------------------------------------------%

% RBF basis function -----------------------------------------------------%
if strmatch('rbf',intype)
    [t1,t2] = meshgrid(time,time);
    in = exp(-(t1-t2).^2./(2*inparam));
end
%-------------------------------------------------------------------------%

% Sigmoidal basis function -----------------------------------------------------%
if strmatch('sigmoid',intype)
    [t1,t2] = meshgrid(time,time);
    in = (1+exp(-(t1-t2)/inparam)).^-1;
end
%-------------------------------------------------------------------------%

% Trangle input ----------------------------------------------------------%
if strmatch('trangle',intype)
    ratio = inparam;
    in = zeros(1,totsamp);
    in(2:ratio:end) = 1;
    for ii = 2:ratio:totsamp
        in(ii:ii+ratio-1)=linspace(1,0,ratio);
    end
end
%-------------------------------------------------------------------------%

% Box input --------------------------------------------------------------%
if strmatch('box',intype)
    ratio = inparam(1);
    in = zeros(1,totsamp);
    in(2:ratio:end) = 1;
    boxwidth = inparam(2);
    for ii = 2:ratio:totsamp
        in(ii:ii+boxwidth-1)=ones(1,boxwidth);
    end
end
%-------------------------------------------------------------------------%

% Exponential decay ------------------------------------------------------%
if strmatch('expdecay',intype)
    ratio = inparam(1);
    in = zeros(1,totsamp);
    in(2:ratio:end) = 1;
    dwidth = inparam(2);
    expdecay = exp(-0.01*dwidth);
    for ii = 2:ratio:totsamp
        in(ii:ii+ratio)=expdecay;
    end
end
%-------------------------------------------------------------------------%

% Bell shape -------------------------------------------------------------%
if strmatch('bell',intype)
    ratio = inparam;
    in = zeros(1,totsamp);
    in(250:ratio:end) = 1;
    width = 0.005:0.005:1;
    bell = exp(-(width-width(end)/2).^2/0.02);
    for ii = 250:ratio:totsamp
        in(ii-ceil(length(width)/2):ii+ceil(length(width)/2))=[0 bell];
    end
end
%-------------------------------------------------------------------------%







