function data = rr2pp(data)
% data = hr2pp(data)
% 
% usage:
%   convert R-R intervals to point processes 
% 
% input:
%   data --- data structure
%
% output
%   data --- processed data
% -------------------------------------
% author: ke yuan 
% email : ky08r@ecs.soton.ac.uk

idx  = 0;                                      
ridx = zeros(1,length(data.rr));
for ii = 1:length(data.rr)
    idx      = idx + data.rr(ii);
    ridx(ii) = idx;
end
data.tottime           = sum(data.rr);
data.totsamp           = fix(data.tottime/data.delta) + 1;
data.time              = 0:data.delta:data.tottime;
ridx                   = fix(ridx/data.delta);
data.observation       = zeros(1,data.totsamp);
data.observation(ridx) = 1;

end