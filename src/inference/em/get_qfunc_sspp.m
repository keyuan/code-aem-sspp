function qfunc = get_qfunc_sspp(y, delta, stats, param, varargin)
% qfunc = get_qfunc_sspp(stats, param, stats, y, delta)
% 
% usage:
%   Compute <log(x,y)>_p(x|y) (aka Q-funciotn) for SSPP model.
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

% Get dim
[totchan,totsamp] = size(y);

% Get parameters
rho     = param.est.rho;
alpha   = param.est.alpha;
sigmasq = param.est.sigmasq;
beta    = param.est.beta;
mu      = param.est.mu;
covinit = param.est.covinit;

% Get sufficient statistics
xsmth          = stats.xsmth;
covsmth        = stats.covsmth;
autoexp        = stats.autoexp;
suminin        = stats.suminin;
sumxcurin      = stats.sumxcurin;
sumxprevin     = stats.sumxprevin;
sumcrossexp    = stats.sumcrossexp;
sumautoexpcur  = stats.sumautoexpcur;
sumautoexpprev = stats.sumautoexpprev;

stadim = size(xsmth,1);

% Compute expected log-complete data likelihood
if (compxinit)
    % Compute <p(x_0)>
    qfunc0 = trace(-covinit\autoexp(:,1)/2) - 2\log(det(covinit));
end

if (compx)
    % Compute <p(x)>
    qfunc1 = trace(-(2*sigmasq)\(sumautoexpcur + rho*sumautoexpprev*rho' + ...
             alpha*suminin*alpha' - rho*sumcrossexp' - sumcrossexp*rho'  - ...
             sumxcurin*alpha' - alpha*sumxcurin' + rho*sumxprevin*alpha' + ...
             alpha*sumxprevin'*rho')) - totsamp/2*log(det(sigmasq)) -      ...
             totsamp*stadim/2*log(2*pi);
end

if (compy)
    % Compute <p(y|x)>
    qfunc2 = 0;
    for nchan = 1:totchan
        bb = beta(:,nchan);
        for nsamp = 1:totsamp
            qfunc2 = qfunc2 + y(nchan,nsamp)*(mu+bb'*xsmth(:,nsamp) +...
                     log(delta)) - exp(mu+bb'*xsmth(:,nsamp) + ...
                     2\bb'*covsmth(:,:,nsamp)*bb)*delta;
        end
    end
end

qfunc  = qfunc0 + qfunc1 +qfunc2;
end



