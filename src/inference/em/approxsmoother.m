function [stats,lb] = approxsmoother(y,in,delta,stats,param,option)
% Purpose:
%   Approximate RTS smoother for SSPP model due to Smith and Brown (2003). 
% Usage:
%   [stats,lb] = approxsmoother(param,y,stats,interm,delta)
% Input:
%   y      --- Binary data
%   stats  --- Estimated sufficient statistics
%   param  --- Estimated parameters
%   interm --- Input related terms
%   delta  --- Time resolution
% Output:
%   stats  --- Estimated sufficient statistics
%   lb     --- Estimated lower bound
%----------------------------------------------
% Author: Ke Yuan 
% Email: ky08r@ecs.soton.ac.uk

% Get option
fltopt = inferget(option,'fltopt');
stadim = inferget(option,'stadim');

% Set default choice for fltopt
if isempty(fltopt), fltopt = 'newton'; end

% Forward pass: Laplace Gaussian filter ----------------------------------%
if stadim == 1
    stats = lgf_sspp_s(y, in, delta, stats, param, fltopt); 
else
    stats = lgf_sspp(y, in, delta, stats, param, fltopt);
end
%-------------------------------------------------------------------------%

% Bacward pass: RTS smoother ---------------------------------------------%
if stadim == 1
    stats = rts_s(stats, param);   
else
    stats = rts(stats, param);   
end
%-------------------------------------------------------------------------%
    
% Get all sufficient statistics related to the joint smoothing density ---%
if stadim == 1
    stats = get_sufstats_s(stats,in);   
else 
    stats = get_sufstats(stats,in);
end
%-------------------------------------------------------------------------%

% Compute current lower bound --------------------------------------------%
if stadim == 1
    qfunc = get_qfunc_sspp_s(y, delta, stats, param);      % Q-Function  
%     ent   = get_ent_sspp_s(stats);                         % Entropy
else
    qfunc = get_qfunc_sspp(y, delta, stats, param);        % Q-Function  
%     ent   = get_ent_sspp(stats);                           % Entropy
end
lb = qfunc;                                                % Lower Bound
%-------------------------------------------------------------------------%


        
