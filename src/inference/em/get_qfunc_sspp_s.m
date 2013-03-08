function qfunc = get_qfunc_sspp_s(y, delta, stats, param, varargin)
% qfunc = get_qfunc_sspp(stats, param, stats, y, delta)
% 
% usage:
%   Compute <log(x,y)>_p(x|y) (aka Q-funciotn) for HESSPP model.
%
% input:
%   stats --- sufficient statistics for joint smoothing density
%   param    --- parameter estimates
%   stats    --- stats estimates
%   y        --- observations
%   delta    --- time resolution 
%   varargin --- optional input chosing likelihood parts  
%
% output:
%   qfunc    --- Q-funciton
%
% -------------------------------------
% author: ke yuan 
% email : ky08r@ecs.soton.ac.uk

if nargin < 5
    compx     = 1;
    compy     = 1;
    compxinit = 0;
elseif (varargin{1} == 1)
    compx     = 1;
    compy     = 0;
    compxinit = 0;
elseif (varargin{1} == 2)
    compx     = 0;
    compy     = 1;
    compxinit = 0;
end

if (compx == 0)
    qfunc1=0;
end
if (compy == 0)
    qfunc2 = 0;
end
if (compxinit == 0)
    qfunc0 = 0;
end

% get dim
[totchan,totsamp] = size(y);

% get parameters
rho     = param.est.rho;
alpha   = param.est.alpha;
sigmasq = param.est.sigmasq;
beta    = param.est.beta';
mu      = param.est.mu;
covinit = param.est.covinit;

% get sufficient statistics
xsmth          = stats.xsmth;
covsmth        = stats.covsmth;
autoexp        = stats.autoexp;
suminin        = stats.suminin;
sumxcurin      = stats.sumxcurin;
sumxprevin     = stats.sumxprevin;
sumcrossexp    = stats.sumcrossexp;
sumautoexpcur  = stats.sumautoexpcur;
sumautoexpprev = stats.sumautoexpprev;

% Compute expected log-complete data likelihood
if (compxinit)
    qfunc0 = -2\autoexp(:,1)/covinit-2\log(covinit);
end

if (compx)
    qfunc1 = -(2*sigmasq)\(sumautoexpcur + rho^2*sumautoexpprev + ...
             alpha*suminin*alpha' - 2*rho*sumcrossexp - ...
             2*sumxcurin*alpha' + 2*rho*sumxprevin*alpha') ...
             -totsamp/2*log(sigmasq) - totsamp/2*log(2*pi);
end

if (compy)
    
    % Vectorize
    xsmth = xsmth(ones(1,totchan),:);
    covsmth = covsmth(ones(1,totchan),:);
    beta = beta(:,ones(1,totsamp-1));
    
    % Compute <p(y|x)>
    qfunc2 = sum(sum(y(:,2:end).*(mu+beta.*xsmth(:,2:end) +...
             log(delta)) - exp(mu+beta.*xsmth(:,2:end) + ...
             beta.^2/2.*covsmth(:,2:end))*delta));
end

qfunc  = qfunc0 + qfunc1 +qfunc2;
end
