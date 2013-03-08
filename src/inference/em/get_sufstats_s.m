function stats = get_sufstats_s(stats,in)
% stats = get_stats_s(stats,in1,in2)
% 
% usage:
%   Compute sufficient statistics for joint smoothing density p(x|y). statss
%   are concerned as scaler.
%
% input:
%   stats --- stats estimates  
%   in    --- input signals
%
% output:
%   stats --- sufficient statistics for joint smoothing density
%
% -------------------------------------
% author: ke yuan 
% email : ky08r@ecs.soton.ac.uk

% get sufficient statistics
xsmth    = stats.xsmth;
autoexp  = stats.autoexp;
crossexp = stats.crossexp;

% compute sumations
suminin        = in(:,2:end)*in(:,2:end)';
sumxcurin      = xsmth(:,2:end)*in(:,2:end)';    
sumxprevin     = xsmth(:,1:end-1)*in(:,2:end)';
sumcrossexp    = sum(crossexp);
sumautoexpcur  = sum(autoexp(:,2:end));    
sumautoexpprev = sum(autoexp(:,1:end-1));    
     
%- align return variables
stats.suminin        = suminin;
stats.sumxcurin      = sumxcurin;
stats.sumxprevin     = sumxprevin;
stats.sumcrossexp    = sumcrossexp;
stats.sumautoexpcur  = sumautoexpcur;
stats.sumautoexpprev = sumautoexpprev;