function in = input_sspp(time,freq,opt)

si = 0;
he = 0;
tr = 0;
box = 0;
expd = 0;
bs   = 0;

if opt == 1; si = 1; end
if opt == 2; he = 1; end
if opt == 3; tr = 1; end
if opt == 4; box = 1; end
if opt == 5; expd = 1; end
if opt == 6; bs = 1; end

totsamp = length(time);

% spike inputs
if (si)
    totin = totsamp/100;
    ratio = ceil(totsamp/totin);  
    in = zeros(1,totsamp);
    in(2:ratio:end) = 1;
end

% harmonic excitation signals
if (he)
    num  = length(freq); 
    in   = zeros(2*num,totsamp);
    for ii = 1:num
        blk = (ii-1)*2+1:ii*2;
        in(blk,:) = [cos(freq(ii)*time);sin(freq(ii)*time)];
    end
end

% trangle input
if (tr)
    totin = totsamp/100;
    ratio = ceil(totsamp/totin);
    in = zeros(1,totsamp);
    in(2:ratio:end) = 1;
    for ii = 2:ratio:totsamp
        in(ii:ii+ratio-1)=linspace(1,0,ratio);
    end
end

% box input
if (box)
    totin = 20;
    ratio = ceil(totsamp/totin);
    in = zeros(1,totsamp);
    in(2:ratio:end) = 1;
    boxwidth = 25;
    for ii = 2:ratio:totsamp
        in(ii:ii+boxwidth-1)=ones(1,boxwidth);
    end
end

% exp decay
if (expd)
    totin = 20;
    ratio = ceil(totsamp/totin);
    in = zeros(1,totsamp);
    in(2:ratio:end) = 1;
    dwidth = 0:ratio;
    expdecay = exp(-0.01*dwidth);
    for ii = 2:ratio:totsamp
        in(ii:ii+ratio)=expdecay;
    end
end

% bell shape
if (bs)
    totin = 20;
    ratio = ceil(totsamp/totin);
    in = zeros(1,totsamp);
    in(250:ratio:end) = 1;
    width = 0.005:0.005:1;
    bell = exp(-(width-width(end)/2).^2/0.02);
    for ii = 250:ratio:totsamp
        in(ii-ceil(length(width)/2):ii+ceil(length(width)/2))=[0 bell];
    end
end







