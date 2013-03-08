function param = mstep_sspp(y, delta,stats, param, fixparam) 
% Purpose:
%   M-step of EM for SSPP model
% Usage:
%   param = mstep_sspp(stats, param, y, fixparam)
% Input:
%   stats --- Sufficient statistics of joint smooth density
%   param --- Estimated parameters
%   y     --- Observations
%   delta --- Time resolution 
% Output:
%   param --- Estimated parameters
% -------------------------------------
% Author: Ke Yuan 
% Email : ky08r@ecs.soton.ac.uk

% Get dimension
[totchan,totsamp] = size(y);

% Get sufficient statistics ----------------------------------------------%
xsmth          = stats.xsmth;
covsmth        = stats.covsmth;
suminin        = stats.suminin;
sumxcurin      = stats.sumxcurin;
sumxprevin     = stats.sumxprevin;
sumcrossexp    = stats.sumcrossexp;
sumautoexpcur  = stats.sumautoexpcur;
sumautoexpprev = stats.sumautoexpprev;
%-------------------------------------------------------------------------%

stadim         = size(xsmth,1);

% Estimate AR(1) coefficient and input weights ---------------------------%
if isempty(strmatch('rho',fixparam)) || isempty(strmatch('alpha',fixparam))
    
    sys     = [sumautoexpprev, sumxprevin;...
               sumxprevin',   suminin];
    out     = [sumcrossexp, sumxcurin];
    unknown = out / sys;

    % Align return variables
    if isempty(strmatch('rho',fixparam))
        rho           = unknown(:,1:stadim);
        param.est.rho = rho;
    end
    if isempty(strmatch('alpha',fixparam))
        alpha = unknown(:,stadim+1:end);
        if stadim == 1, alpha = alpha'; end
        param.est.alpha = alpha;
    end
end
%-------------------------------------------------------------------------%

% Estimate AR(1) noise variance ------------------------------------------%
if isempty(strmatch('sigma',fixparam)) 
    rho = param.est.rho; 
    alpha = param.est.alpha;
    if ~isempty(strmatch('alpha',fixparam)) ||                        ... 
            ~isempty(strmatch('rho',fixparam))   
        statsterms = sumautoexpcur + rho*sumautoexpprev*rho' +        ...
                     alpha*suminin*alpha' - rho*sumcrossexp' -        ... 
                     sumcrossexp*rho' - sumxcurin*alpha' -            ...
                     alpha*sumxcurin' + rho*sumxprevin*alpha' +       ...
                     alpha*sumxprevin'*rho';
    else
        statsterms = sumautoexpcur - rho*sumcrossexp' - alpha*sumxcurin';
    end
    sigmasq = totsamp\diag(diag(statsterms));
    
    % Align return variables
    param.est.sigmasq = sigmasq;
end
%-------------------------------------------------------------------------%

% Estimate state weight --------------------------------------------------%
if isempty(strmatch('beta',fixparam))
    if stadim == 1
        [mu,beta] = newton_s(param,y,xsmth,covsmth,delta); 
    else
        [mu,beta] = newton(param,y,xsmth,covsmth,delta);
    end
    if isempty(strmatch('mu',fixparam))
        param.est.mu = mean(mu);
    end
    param.est.beta = beta;
end
%-------------------------------------------------------------------------%

% Estimate background firing rate ----------------------------------------%
if isempty(strmatch('mu',fixparam)) && ~isempty(strmatch('beta',fixparam))
    mu = zeros(1,totchan);
    if stadim == 1
        for nchan = 1 : totchan
            mu(nchan) = log(sum(y(nchan,:))) -                       ...
                        log(sum(exp(param.est.beta(nchan)*xsmth +    ...
                        2\param.est.beta(nchan)^2*covsmth)*delta));
        end
    else
        for nchan = 1:totchan
            t1 = 0; bb = param.est.beta(:,nchan);
            for nsamp = 1:totsamp
                x = xsmth(:,nsamp); cov = covsmth(:,:,nsamp); 
                t1 = t1 + exp(bb'*x + 2\bb'*cov*bb)*delta;
            end
            mu(nchan) = log(sum(y(nchan,:))) - log(t1); 
        end
    end
    
    mu = mean(mu);

    % Align return variables
    param.est.mu = mu;
end
%-------------------------------------------------------------------------%

% Update initial stats and its covariance --------------------------------%
if stadim == 1
    xinit   = param.est.rho*stats.xsmth(:,2); 
    covinit = param.est.sigmasq/(1-param.est.rho^2);
else 
    xinit = xsmth(:,1);
    covinit = covsmth(:,:,1);
end

% Align return variables
param.est.xinit   = xinit;
param.est.covinit = covinit;
%-------------------------------------------------------------------------%

end

% SUBROUTINES ============================================================%

function [mu,beta] = newton_s(param,y,xsmth,covsmth,delta)
[totchan,~] = size(y);
mu = repmat(param.est.mu,totchan,1);
beta = param.est.beta;
totnewton = 6;

for nchan = 1:totchan
    theta = [mu(nchan);beta(nchan)];    
    for nnewton = 1:totnewton
        mm = mu(nchan); bb = beta(nchan);
        t1 = exp(mm + bb*xsmth +  ...
                    2\bb^2*covsmth)*delta;
        t2 = xsmth+covsmth*bb;
        g1 = sum(y(nchan,:) - t1 );
        g2 = sum(y(nchan,:).*xsmth-t2.*t1);
        h1 = -sum(t1);
        h3 = -sum(t2.*t1);
        h2 = h3;
        h4 = -sum(covsmth.*t1+t2.^2.*t1);
        g  = [g1;g2]; H = [h1,h2;h3,h4];
        theta = theta -H\g;
        mu(nchan) = theta(1); beta(nchan) = theta(2);
    end
end
end

function [mu,beta] = newton(param,y,xsmth,covsmth,delta)
[totchan,totsamp] = size(y);
mu = repmat(param.est.mu,totchan,1);
beta = param.est.beta;
totnewton = 6;

for nchan = 1:totchan
    theta = [mu(nchan);beta(:,nchan)];
    for nnewton = 1:totnewton
        g1 = 0; g2 = 0; h1 = 0; h3 = 0; h4 = 0;
        mm = mu(nchan); bb = beta(:,nchan);
        for nsamp = 1:totsamp
            x = xsmth(:,nsamp);
            cov = covsmth(:,:,nsamp);
            t1 = exp(mm + bb'*x + 2\bb'*cov*bb)*delta;
            t2 = x+cov*bb;
            g1 = g1 + y(nchan,nsamp) - t1;
            g2 = g2 + y(nchan,nsamp)*x - t2*t1;
            h1 = h1 - t1;
            h3 = h3 - t2*t1;
            h2 = h3';
            h4 = h4 - (cov*t1 + t2*t1*t2');
        end
        g  = [g1;g2]; H = [h1,h2;h3,h4];
        theta = theta - H\g;
        mu(nchan) = theta(1); beta(:,nchan) = theta(2:end);
    end
end
end
