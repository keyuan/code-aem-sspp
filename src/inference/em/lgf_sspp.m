function stats = lgf_sspp(y, in, delta, stats, param, fltopt)
% Purpose: 
%   Lapace Gaussian filter (LGF) for state-space model with point process
%   observations (SSPP). States are concerned as vecoters.
% Usage:
%   stats = lgf_sspp_s(param, y, stats, interm, delta, fltopt)
% Input:
%   param  --- Parameters of SSPP model
%   y      --- Observations  
%   interm --- Externel stimuli terms (alpha*I_k)
%   fltopt --- Options for chosing a optimiser to obtain MAP solution 
%                of filtering density 
% Output:
%   xpost   --- Filtering density mode
%   covpost --- Filtering density covariance 
% ------------------------------------------
% Author: Ke Yuan 
% Email : ky08r@ecs.soton.ac.uk

% Preallocate variables for speed
xpred    = stats.xpred;                  
xpost    = stats.xpost;
covpred  = stats.covpred;          
covpost  = stats.covpost;

% Get dimensions
[totchan,totsamp] = size(y);

% Get parameters
rho     = param.est.rho;
sigmasq = param.est.sigmasq;
beta    = param.est.beta;
mu      = param.est.mu;
xinit   = param.est.xinit;
covinit = param.est.covinit;

% Set initial conditions
xpost(:,1)   = xinit;
covpost(:,:,1) = covinit;

% Input terms
interm = param.est.alpha*in;
idm    = eye(size(xinit,1));

% FORWARD LOOP ===========================================================%
for nsamp = 2 : totsamp 
    
    %= One-step prediction
    xpred(:,nsamp) = rho*xpost(:,nsamp-1) + interm(:,nsamp);
    
    %= One-step prediction covariance
    covpred(:,:,nsamp) = rho*covpost(:,:,nsamp-1)*rho' + sigmasq;                            
    
    % Compute posterior mode via Newton's method -------------------------% 
    if strcmp('newton',fltopt), totnewton = 10;
        xpost(:,nsamp) = newton(y(:,nsamp), xpred(:,nsamp), ...
            covpred(:,:,nsamp), mu, beta, delta, totnewton);
    end
    %= -------------------------------------------------------------------%
    
    % Compute posterior mode via fix point method ------------------------%
    if strcmp('fixpt',fltopt), totfixpt = 6;
         xpost(:,nsamp) = fixpt(y(:,nsamp), xpred(:,nsamp), ...
             covpred(:,:,nsamp), mu, beta, delta, totfixpt);
    end
    %= -------------------------------------------------------------------%
    
    % Computing posterior mode via replace rhs's x_k|k with x_k|k-1 ------%
    if strcmp('simple',fltopt)
        xpost(:,nsamp) = simple(y(:,nsamp), xpred(:,nsamp), ...
            covpred(:,:,nsamp), mu, beta, delta);
    end 
    %---------------------------------------------------------------------%
    
    t2 = 0;
    for nchan = 1:totchan
        bb = beta(:,nchan);
        t1 = exp(mu+bb'*xpost(:,nsamp))*delta;
        t2 = t2 - bb*t1*bb';
    end
    covpost(:,:,nsamp) = (covpred(:,:,nsamp)\idm - t2)\idm;      
end
%=========================================================================%

% Align return variables
stats.xpred   = xpred;                  
stats.xpost   = xpost;
stats.covpred = covpred;          
stats.covpost = covpost;
end


% SUBROUTINES ============================================================%

function  xpost = newton(y,xpred,covpred,mu,beta,delta,totnewton)
% Compute posterior mode via Newton's method -----------------------------% 
totchan = size(y,1); px = xpred;
for nnewton = 1:totnewton
    t1 = 0; t2 = 0;
    for nchan = 1:totchan
        bb   = beta(:,nchan);
        cifd = exp(mu+bb'*px)*delta;
        t1   = t1 + bb*(y(nchan)-cifd);    
        t2   = t2 - bb*cifd*bb';    
    end
    fx   = xpred+covpred*t1-px;
    dfx  = covpred*t2 - eye(size(xpred,1));
    px   = px - dfx\fx;
end
xpost = px;
end
%= -----------------------------------------------------------------------%

function xpost = fixpt(y,xpred,covpred,mu,beta,delta,totfixpt)
% Compute posterior mode via fix point method ----------------------------%
totchan = size(y,1); px = xpred;
for nfixpt = 1:totfixpt
    t1 = 0;
    for nchan = 1:totchan
        bb = beta(:,nchan);
        t1 = t1 + bb*(y(nchan)-delta*exp(mu+bb'*px));
    end
    px = xpred + covpred*t1;
    
end
xpost = px;
end
%-------------------------------------------------------------------------%

function xpost = simple(y,xpred,covpred,mu,beta,delta)
% Computing posterior mode via replace rhs's x_k|k with x_k|k-1 ----------%
totchan = size(y,1); t1 = 0;
for nchan = 1:totchan
    bb = beta(:,nchan);    
    t1 = t1 + bb*(y(nchan)-delta*exp(mu+bb'*xpred));   
end
xpost = xpred + covpred*t1;

end
%-------------------------------------------------------------------------%

%=========================================================================%
