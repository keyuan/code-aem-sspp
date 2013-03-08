function synthdata = synthdatapp(dim, param, option)
% syndata = he_sspp_syntheticdata(param, dim)
% 
% Purpose:
%   Generate synthetic data from point process models. 
% 
% input:
%   dim    --- dimension settings
%   param  --- parameter settings
%   option --- option from INFERSET
%
% output:
%   synthdata --- resulting data
% -------------------------------------
% author: ke yuan 
% email : ky08r@ecs.soton.ac.uk

% Get option 
model = inferget(option,'model','default','fast');
intype = inferget(option,'intype','default','fast');
ciftype = inferget(option,'cif','default','fast');

% Get dimensions
delta   = dim.delta;                           % Time resolution                               
totchan = dim.totchan;                         % Number of channels
stadim  = dim.stadim;                          % Dimension of state
tottime = dim.tottime;                         % Total observation time
time    = 0:delta:tottime;                     % Time axis
totsamp = tottime/delta + 1;                   % Total sample points
tothist = length(param.true.gamma);
inparam = dim.inparam;

% Get true parameter 
mu     = param.true.mu;                        % Background firing rate
if length(mu) == 1
    mu = repmat(mu,totchan,1);
end
gamma = param.true.gamma;                     % Weight of history 
gamma = repmat(gamma,1,totchan);

% preallocate variables
y   = zeros(totchan,totsamp);
cif = zeros(totchan,totsamp);

% get excitation signals
if strcmp(intype,'none')
    in = zeros(1,totsamp);
else
in  = input_sspp(time,intype,inparam);
end
% initial state
xinit = param.true.xinit;

% Generate hidden states -------------------------------------------------% 
if strcmp(model,'sspp')
    rho     = param.true.rho;                       % AR coefficient
    alpha   = param.true.alpha;
    sigmasq = param.true.sigmasq;                 % State noise variance
    covinit = param.true.covinit;
    beta    = param.true.beta;
    sigma   = sqrt(sigmasq);
    x       = zeros(stadim,totsamp);
    x(:,1)  = xinit + sqrt(covinit)*randn(size(xinit));
    % Generate states
    for nsamp = 2:totsamp
       if stadim == 1
           alpha = alpha';
       end
       x(:,nsamp) = rho*x(:,nsamp-1) + alpha*in(:,nsamp) + ...
                    sigma*randn(size(xinit));         
    end    
end
%-------------------------------------------------------------------------%


% Vectorize parameter vector for PPGLM
if strcmp(model,'ppglm')
    x = repmat(xinit,1,totsamp);
end
%-------------------------------------------------------------------------%

% Generate observaitons via PPGLM ----------------------------------------%
if strcmp(model,'ppglm')
    if isempty(gamma)
        for nsamp = 1:totsamp
            % CIF (Rate Funciton) input
            cifargin = mu+in(:,nsamp)'*x(:,nsamp);
            % Assign pulses
            [y(:,nsamp), cif(:,nsamp)] = inonlygenerator...
                (ciftype,cifargin,delta);
        end
    else % Use history information
        % Initial CIF and pulses ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        cifargin = mu+sum(in(:,1:tothist).*x(:,1:tothist),1);    
        [y(:,1:tothist), cif(1:tothist)] = inonlygenerator...
            (ciftype,cifargin,delta);
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%    
        for nsamp = tothist+1:totsamp
            % CIF (Rate Funciton)
            cifargin = mu+in(:,nsamp)'*x(:,nsamp);
            yhist = y(:,nsamp-1:-1:nsamp-tothist)';
            % Assign pulses
            [y(:,nsamp), cif(:,nsamp)] = inhistgenerator...
               (ciftype,yhist,gamma,cifargin,delta);
        end
        
    end
end
%-------------------------------------------------------------------------%

% Generate observaitons via SSPP -----------------------------------------%
if strcmp(model,'sspp')
    for nsamp = 1:totsamp
        for nchan = 1:totchan
            % CIF (Rate Funciton) input
            cifargin = mu(nchan)+beta(:,nchan)'*x(:,nsamp);
            % Assign pulses
            [y(nchan,nsamp), cif(nchan,nsamp)] = inonlygenerator  ...
                (ciftype,cifargin,delta);
        end
    end
end
%-------------------------------------------------------------------------%

% Align return variables
synthdata = struct('type', 'synthetic data from piont process');
synthdata.x       = x;
synthdata.in      = in;
synthdata.cif     = cif;
synthdata.y       = y;
synthdata.time    = time;
synthdata.delta   = delta;
synthdata.inparam = inparam;
end

function [y, cif] = inonlygenerator(ciftype,cifargin,delta)
% Generate pulses without history information
cif = feval(ciftype,cifargin);
if strcmp(ciftype,'sigmoid'),cifd = cif; 
else cifd = cif*delta; end
y   = (rand(size(cif)) < cifd);
end

function [y,cif] = inhistgenerator(ciftype,yhist,gamma,cifargin,delta)
% Generate pulses with input and history information
cif = feval(ciftype,cifargin + gamma'*yhist);
if strcmp(ciftype,'sigmoid'),cifd = cif; 
else cifd = cif*delta; end
y   = (rand(size(cif)) < cifd);
end

function f = sigmoid(z)
f = (exp(-z)+1).^-1;
end


