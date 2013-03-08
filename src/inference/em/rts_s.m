function stats = rts_s(stats, param)
% stats = rts_s(stats)
%
% usage:
%   RTS (fixed interval) smoother for scaler stats systems
%
% input:
%   stats --- stats estimates from previous step
%   param --- parameters 
%
% output:
%   stats --- smoothed stats estimates
%
% -------------------------------------
% author: ke yuan 
% email : ky08r@ecs.soton.ac.uk

% get stats estimates
xpred    = stats.xpred;
covpred  = stats.covpred;
xpost    = stats.xpost;
covpost  = stats.covpost;
xsmth    = stats.xsmth;
covsmth  = stats.covsmth;
crosscov = stats.crosscov;
autoexp  = stats.autoexp;
crossexp = stats.crossexp;

% get dimensions
totsamp = size(xpred,2);

% get parameter
rho = param.est.rho;

% backward pass: RTS smoothing -------------------------------------------%

%- initialize variables 
xsmth(:,end)   = xpost(:,end);                                                                    
covsmth(:,end) = covpost(:,end);
autoexp(:,end) = covsmth(:,end) + xsmth(:,end)*xsmth(:,end)';

%- backward loop
for nsamp = totsamp-1:-1:1
    jk                = covpost(:,nsamp)*rho'/covpred(:,nsamp+1);  
    xsmth(:,nsamp)    = xpost(:,nsamp) + ...
                        jk*(xsmth(:,nsamp+1)-xpred(:,nsamp+1));
    covsmth(:,nsamp)  = covpost(:,nsamp) + ...
                        jk*(covsmth(:,nsamp+1)-covpred(:,nsamp+1))*jk';            
    autoexp(:,nsamp)  = covsmth(:,nsamp) + xsmth(:,nsamp)*xsmth(:,nsamp)';
    crosscov(:,nsamp) = jk*covsmth(:,nsamp+1);
    crossexp(:,nsamp) = crosscov(:,nsamp) + xsmth(:,nsamp+1)*xsmth(:,nsamp)';
end
%-------------------------------------------------------------------------%

% align return variables
stats.xsmth    = xsmth;
stats.covsmth  = covsmth;
stats.crosscov = crosscov;
stats.autoexp  = autoexp;
stats.crossexp = crossexp;

end