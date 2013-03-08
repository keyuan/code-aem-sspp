function data = pp2rr(data)
% data = hr2pp(data)
% 
% usage:
%   convert point processes to R-R intervals
% 
% input:
%   data --- data structure
%
% output
%   data --- processed data
% -------------------------------------
% author: ke yuan 
% email : ky08r@ecs.soton.ac.uk

y       = data.observation;
idx     = find(y==1);
drr     = idx(2:end)-idx(1:end-1);
rr      = drr*data.delta;
data.rr = rr;