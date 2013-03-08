function stats = lgf_sspp_s(y, in, delta, stats, param, fltopt)
% Purpose: 
%   Lapace Gaussian filter (LGF) for state-space model with point process
%   observations (SSPP). States are concerned as scaler.
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
totsamp  = size(xpred,2);

% Get parameters
rho     = param.est.rho;
sigmasq = param.est.sigmasq;
beta    = param.est.beta;
mu      = param.est.mu;
xinit   = param.est.xinit;
covinit = param.est.covinit;

% Set initial conditions
xpost(:,1)   = xinit;
covpost(:,1) = covinit;

% Input terms
interm = param.est.alpha*in;

% FORWARD LOOP ===========================================================%
for nsamp = 2 : totsamp 
    
    %= One-step prediction
    xpred(:,nsamp) = rho*xpost(:,nsamp-1) + interm(:,nsamp);
    
    %= One-step prediction covariance
    covpred(:,nsamp) = rho*covpost(:,nsamp-1)*rho' + sigmasq;                            
    
    % Compute posterior mode via Newton's method -------------------------% 
    if strcmp('newton',fltopt), totnewton = 5;
        xpost(:,nsamp) = newton(y(:,nsamp), xpred(:,nsamp), ...
            covpred(:,nsamp), mu, beta, delta, totnewton);
    end
    %= -------------------------------------------------------------------%
    
    % Compute posterior mode via fix point method ------------------------%
    if strcmp('fixpt',fltopt), totfixpt = 5;
         xpost(:,nsamp) = fixpt(y(:,nsamp), xpred(:,nsamp), ...
             covpred(:,nsamp), mu, beta, delta, totfixpt);
    end
    %= -------------------------------------------------------------------%
    
    % Computing posterior mode via replace rhs's x_k|k with x_k|k-1 ------%
    if strcmp('simple',fltopt)
        xpost(:,nsamp) = simple(y(:,nsamp), xpred(:,nsamp), ...
            covpred(:,nsamp), mu, beta, delta);
    end 
    %---------------------------------------------------------------------%
    
    hessianlikhd     = -sum(beta.^2.*exp(mu + beta*xpost(:,nsamp))*delta);
    covpost(:,nsamp) = (covpred(:,nsamp)\1 - hessianlikhd)\1;      
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
px = xpred;
for nnewton = 1:totnewton
    gradlikhd    = sum(beta.*(y'-delta*exp(mu+beta*px))); 
    fx           = xpred+covpred*gradlikhd-px;
    hessianlikhd = -sum(beta.^2.*exp(mu+beta*px)*delta);
    dfx          = covpred*hessianlikhd-1;
    px           = px - fx/dfx;   
end
xpost = px;
end
%= -----------------------------------------------------------------------%

function xpost = fixpt(y,xpred,covpred,mu,beta,delta,totfixpt)
% Compute posterior mode via fix point method ----------------------------%
px = xpred;
for nfixpt = 1:totfixpt
    gradlikhd = sum(beta.*(y'-delta*exp(mu+beta*px)));
    px        = xpred + covpred*gradlikhd;
end
xpost = px;
end
%-------------------------------------------------------------------------%

function xpost = simple(y,xpred,covpred,mu,beta,delta)
% Computing posterior mode via replace rhs's x_k|k with x_k|k-1 ----------%
gradlikhd = sum(beta.*(y'-delta*exp(mu+beta*xpred)));
xpost = xpred + covpred*gradlikhd;
end
%-------------------------------------------------------------------------%

%=========================================================================%
