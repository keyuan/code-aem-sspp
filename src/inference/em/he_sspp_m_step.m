function param = he_sspp_m_step(stats, param, y, delta, varargin)
% param = he_sspp_m_step(stats, param, stats, y, delta)
%
% usage:
%   M-step of EM for SSPP model
% 
% input:
%   stats --- sufficient statistics of joint smooth density
%   param    --- estimated parameters
%   stats    --- stats estimates
%   y        --- observations
%   delta    --- time resolution 
%
% output:
%   param --- estimated parameters
%
% -------------------------------------
% author: ke yuan 
% email : ky08r@ecs.soton.ac.uk

if (nargin < 7)
    compx     = 1;
    compy     = 1;
    compxinit = 1;
elseif (varargin{1} == 1)
    compx     = 1;
    compy     = 0;
    compxinit = 1;
elseif (varargin{1} == 2)
    compx     = 0;
    compy     = 1;
    compxinit = 0;
elseif (varargin{1} == 3)
    compx     = 0;
    compy     = 1;
    compxinit = 1;
end

% get fixed parameter
beta = param.est.beta;

% get dimension
[totchan,totsamp] = size(y);

% get sufficient statistics
xsmth          = stats.xsmth;
covsmth        = stats.covsmth;
suminin        = stats.suminin;
sumxcurin      = stats.sumxcurin;
sumxprevin     = stats.sumxprevin;
sumcrossexp    = stats.sumcrossexp;
sumautoexpcur  = stats.sumautoexpcur;
sumautoexpprev = stats.sumautoexpprev;

if (compx)
    
    % estimate AR(1) coefficient and two input weights
    sys     = [sumautoexpprev, sumxprevin;...
               sumxprevin',   suminin];
    out     = [sumcrossexp; sumxcurin'];
    unknown = sys \ out;
    rho     = unknown(1);
    alpha   = unknown(2:end)';
    
    % align return variables
    param.est.rho   = rho;
    param.est.alpha = alpha;

    % estimate AR(1) noise variance
    fixsigma = 1;
    if (~fixsigma) 
        statsterms = sumautoexpcur + rho^2*sumautoexpprev + ...
                     alpha*suminin*alpha' - 2*rho*sumcrossexp - ...
                     2*sumxcurin*alpha' + 2*rho*sumxprevin*alpha' +...
                     stats.autoexp(:,1)*(1-rho^2);
        sigmasq    = (totsamp-1)\statsterms;
        
        % align return variables
        param.est.sigmasq = sigmasq;
    end
end

if (compy)
    
    % estimate background firing rate
    mu = zeros(1,totchan);
    for nchan = 1 : totchan
        mu(nchan) = log(sum(y(nchan,:))) -...
                    log(sum(exp(beta(nchan,:)*xsmth +...
                    2\beta(nchan,:)^2*covsmth)*delta));
    end
    mu = mean(mu);
    
    % align return variables
    param.est.mu = mu;
end

if (compxinit)
    
    % update initial stats and its covariance
    xinit   = rho*stats.xsmth(:,2); 
    covinit = param.est.sigmasq/(1-param.est.rho^2);
    
    % align return variables
    param.est.xinit   = xinit;
    param.est.covinit = covinit;
end

end
