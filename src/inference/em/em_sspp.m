function [param,stats,lbsave,nem] = em_sspp(data,param,option)
% Purpose:
%   Expectation-Maximization (EM) algorithm for SSPP model
% Usage:
%   [param,stats,lbsave] = em(data,option,param)
% Input:
%   data   --- Observed data
%   option --- Option from INFERSET
%   param  --- Initial conditions
% Output:
%   param  --- Estimated parameters
%   stats  --- Estimated statistics
%   lbsave --- Lower bound trajectry  
%   nem    --- Stoping step
% -------------------------------------
% Author: Ke Yuan 
% Email : ky08r@ecs.soton.ac.uk

% INFR STATES AND PARAMETERS VIA EM ======================================%

% Get options ------------------------------------------------------------%
estep    = inferget(option,'estep');
tol      = inferget(option,'tolfun');
intype   = inferget(option,'intype');
stadim   = inferget(option,'stadim');
totem    = inferget(option,'maxiter');
display  = inferget(option,'display');
fixparam = inferget(option,'fixparam');
%-------------------------------------------------------------------------%

% Get data ---------------------------------------------------------------%
y      = data.y;
time   = data.time;
delta  = data.delta;
if ~isempty(data.in)
    in = data.in;
else
    inparam = data.inparam;
    in   = input_sspp(time,intype,inparam);
end
%-------------------------------------------------------------------------%

% Get dimensions ---------------------------------------------------------%
totsamp = size(y,2);
lbsave  = zeros(1,totem);
%-------------------------------------------------------------------------%

% Preallocate saved parameter estimaiton result --------------------------%
if stadim == 1;  
    param.save.rho    = zeros(stadim,totem);
    param.save.alpha  = zeros(length(param.est.alpha),totem);
    param.save.mu     = zeros(1,totem);
    param.save.beta   = zeros(size(y,1),totem);
else
    param.save.rho    = zeros(stadim*stadim,totem);
    param.save.alpha  = zeros(stadim*length(param.est.alpha),totem);
    param.save.mu     = zeros(1,totem);
    param.save.beta   = zeros(stadim*size(y,1),totem);
end
%-------------------------------------------------------------------------%

% Align initial guesses --------------------------------------------------%
param.save.rho(:,1)   = param.est.rho(:);
param.save.alpha(:,1) = param.est.alpha(:);
param.save.mu(:,1)    = param.est.mu(:);
param.save.beta(:,1)  = param.est.beta(:);
%-------------------------------------------------------------------------%

% Preallocate E-step variables -------------------------------------------%
stats          = struct('type','inferred statistics');
if stadim == 1
    stats.xpred    = zeros(stadim,totsamp);
    stats.covpred  = zeros(stadim,totsamp);
    stats.xpost    = zeros(stadim,totsamp);
    stats.covpost  = zeros(stadim,totsamp);
    stats.xsmth    = zeros(stadim,totsamp);
    stats.covsmth  = zeros(stadim,totsamp);
    stats.autoexp  = zeros(stadim,totsamp);
    stats.crosscov = zeros(stadim,totsamp-1);
    stats.crossexp = zeros(stadim,totsamp-1);
else
    stats.xpred    = zeros(stadim,totsamp);
    stats.covpred  = zeros(stadim,stadim,totsamp);
    stats.xpost    = zeros(stadim,totsamp);
    stats.covpost  = zeros(stadim,stadim,totsamp);
    stats.xsmth    = zeros(stadim,totsamp);
    stats.covsmth  = zeros(stadim,stadim,totsamp);
    stats.autoexp  = zeros(stadim,stadim,totsamp);
    stats.crosscov = zeros(stadim,stadim,totsamp-1);
    stats.crossexp = zeros(stadim,stadim,totsamp-1);
end
%-------------------------------------------------------------------------%

% MAIN LOOP FOR EM =======================================================%
for nem = 2:totem    
    
    % E-Step: Infer smooth density and related statistics ----------------%   
    [stats,lb] = feval(estep, y, in, delta, stats, param, option); 
    lbsave(:,nem) = lb;                      % Save lower bound trajectry
    %---------------------------------------------------------------------%
    
    % M-Step: MLE for parameters -----------------------------------------%
    param = mstep_sspp(y, delta,stats, param, fixparam);       
    %---------------------------------------------------------------------%
    
    % Save learning trajectry --------------------------------------------%
    param.save.rho(:,nem)     = param.est.rho(:);
    param.save.alpha(:,nem)   = param.est.alpha(:);
    param.save.sigmasq(:,nem) = param.est.sigmasq(:);
    param.save.mu(:,nem)      = param.est.mu(:);
    param.save.beta(:,nem)    = param.est.beta(:);
    %---------------------------------------------------------------------%
        
    % Stopping condition -------------------------------------------------%      
    reldiff = abs(lb-lbsave(nem-1))/abs(lbsave(nem-1));
    if (nem > 3)
        if ( reldiff < tol ) || ~isfinite(lb)        
            if strcmp(display,'final')
                fprintf([' EM = %i, Lower Bound = %f, Relative Change '...
                    '= %f\n'], nem, lb, reldiff);
                disp('Estimated Parameters:')
                disp(param.est)
            end
        break
        end
    end
    %---------------------------------------------------------------------%
    
    % Display setting ----------------------------------------------------%
    if strcmp(display,'iter')      
        fprintf(' EM = %i, Lower Bound = %f, Relative Change = %f\r', ...
            nem, lb, reldiff)                   
    end
    %---------------------------------------------------------------------%
    
end
%=========================================================================%
fprintf('\n')
end
