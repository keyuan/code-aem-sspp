function [param,state,lbsave,nem] = he_sspp_em...
                                    (data,stadim,totem,tol,param,varargin)
% [param,state,lbsave] = he_sspp_em(data,stadim,totem,tol,param,varargin)
%
% usage:
%   Expectation-Maximization (EM) algorithm for HESSPP model
%
% input:
%   data     --- observed data
%   stadim   --- dimension of state
%   totem    --- maximum iterations
%   tol      --- convergenc tollrence
%   param    --- initial conditions
%   varargin --- input type
%
% output:
%   param  --- estimated parameters
%   state  --- estimated states
%   lbsave --- lower bound trajectry  
%   nem    --- stoping step
%
% -------------------------------------
% author: ke yuan 
% email : ky08r@ecs.soton.ac.uk



% Infer x and parameters via EM ==========================================%

% get data
y     = data.observation;
time  = data.time;
delta = data.delta;

if (nargin < 6)
    in = data.in;
else
    freq = param.est.freq;
    in   = input_sspp(time,freq,varargin{1});
end

% get dimensions
totsamp = size(y,2);
lbsave  = zeros(1,totem);
const   = totsamp/2;

%= preallocate saved parameter estimaiton result
param.save.rho    = zeros(stadim,totem);
param.save.alpha  = zeros(length(param.est.alpha),totem);
param.save.mu     = zeros(1,totem);

%= align initial guesses 
param.save.rho(:,1)     = param.est.rho;
param.save.alpha(:,1)   = param.est.alpha';
param.save.mu(:,1)      = param.est.mu;

%= preallocate E-step variables
state          = struct('type','state estimates');
state.xpred    = zeros(stadim,totsamp);
state.covpred  = zeros(stadim,totsamp);
state.xpost    = zeros(stadim,totsamp);
state.covpost  = zeros(stadim,totsamp);
state.xsmth    = zeros(stadim,totsamp);
state.covsmth  = zeros(stadim,totsamp);
state.autoexp  = zeros(stadim,totsamp);
state.crosscov = zeros(stadim,totsamp-1);
state.crossexp = zeros(stadim,totsamp-1);

%= initial input terms
interm = param.est.alpha*in;

%= main loop for EM 
for nem = 2:totem
    
    %= E-step ------------------------------------------------------------%
   
    %=- forward pass: laplace gaussian filter
    state = lgf_sspp_s(param, y, state, interm, delta);  
    
    %=- bacward pass: rts smoother
    state = rts_s(state, param);   
    
    %=- get sufficient statistics for the joint smoothing density
    sufstats = get_sufstats_s(state,in);   
    
    %=- compute current lower bound
    qfunc = get_qfunc_sspp(sufstats,param,state,y,delta);    % Q-function  
    ent   = get_ent_sspp_s(state);                           % entropy     
    lb    = qfunc + ent + const...                           % lower bound
            + param.est.sigmasq/(1-param.est.rho^2);
    %=- save lower bound trajectry
    lbsave(:,nem) = lb; 
    %=- ------------------------------------------------------------------%
    
    %= M-Step ------------------------------------------------------------%
    
    %=- parameter estimates
    param = he_sspp_m_step(sufstats, param, state, y, delta); 
    
    %=- update input terms
    interm = param.est.alpha*in;
    %=- ------------------------------------------------------------------%
    
    %= save the estimation trajectry
    param.save.rho(:,nem)     = param.est.rho;
    param.save.alpha(:,nem)   = param.est.alpha';
    param.save.sigmasq(:,nem) = param.est.sigmasq;
    param.save.mu(:,nem)      = param.est.mu;
    
    %= show progress
    if (mod(nem, 10) == 0)        
        disp('-------------------------------------------------------')
        disp(['EM iteration (maximum ' num2str(totem) '):'])
        disp(nem)
        disp('estimated lower bound: ')
        disp(lb)
        disp('lower bound difference:')
        if (nem>10)
            past = 10;
        else past = 8;
        end
        mdf = mean(lbsave(:,nem-past+1:nem)-lbsave(:,nem-past:nem-1));
        disp(mdf)
        disp('estimated parameters:')
        disp(param.est)    
    end
    
    %= stopping condition        
    if (nem > 3)
        diff = abs(lb-lbsave(nem-1))/abs(lbsave(nem-1));
        if ( diff < tol ) || ~isfinite(lb)        
        disp('-------------------------------------------------------')        
        disp([ 'reached convergence at the ' num2str(nem-1) 'th EM step,' ])        
        disp([ 'with relative change precision at ' num2str(tol) '.'])      
        disp('-------------------------------------------------------')       
        break
        end
    end
end
%= =======================================================================%

end
