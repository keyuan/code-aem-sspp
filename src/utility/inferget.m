function o = inferget(options,name,default,flag)
%inferget Get infer OPTIONS parameters.
%   VAL = inferget(OPTIONS,'NAME') extracts the value of the named parameter
%   from inferization options structure OPTIONS, returning an empty matrix if
%   the parameter value is not specified in OPTIONS.  It is sufficient to
%   type only the leading characters that uniquely identify the
%   parameter.  Case is ignored for parameter names.  [] is a valid OPTIONS
%   argument.
%
%   VAL = inferget(OPTIONS,'NAME',DEFAULT) extracts the named parameter as
%   above, but returns DEFAULT if the named parameter is not specified (is [])
%   in OPTIONS.  For example
%
%     val = inferget(opts,'TolX',1e-4);
%
%   returns val = 1e-4 if the TolX parameter is not specified in opts.
%
%   See also INFERSET.

%   Copyright 1984-2008 The MathWorks, Inc.
%   $Revision: 1.20.4.14 $  $Date: 2008/12/01 07:17:28 $

% undocumented usage for fast access with no error checking
if (nargin == 4) && isequal(flag,'fast')
    o = infergetfast(options,name,default);
    return
end

if nargin < 2
    error('inferpp:inferget:NotEnoughInputs', 'Not enough input arguments.');
end
if nargin < 3
    default = [];
end

if ~isempty(options) && ~isa(options,'struct')
    error('inferpp:inferget:Arg1NotStruct',...
        'First argument must be an options structure created with INFERSET.');
end

if isempty(options)
    o = default;
    return;
end

allfields = {'model'; 'likhd'; 'cif'; 'method'; ...
             'maxiter'; 'display'; 'tolfun';'tolarg';'stadim';...
             'fixparam';'intype';'estep';'fltopt'};


Names = allfields;

name = deblank(name(:)'); % force this to be a row vector
j = find(strncmpi(name,Names,length(name)));
if isempty(j)               % if no matches
    error('inferpp:inferget:InvalidPropName',...
        ['Unrecognized option name ''%s''.  ' ...
        'See INFERSET for possibilities.'], name);
elseif length(j) > 1            % if more than one match
    % Check for any exact matches (in case any names are subsets of others)
    k = find(strcmpi(name,Names));
    if length(k) == 1
        j = k;
    else
        msg = sprintf('Ambiguous option name ''%s'' ', name);
        msg = [msg '(' Names{j(1),:}];
        for k = j(2:length(j))'
            msg = [msg ', ' Names{k,:}];
        end
        msg = [msg, '.)'];
        error('inferpp:inferget:AmbiguousPropName', msg);
    end
end

if any(strcmp(Names,Names{j,:}))
    o = options.(Names{j,:});
    if isempty(o)
        o = default;
    end
else
    o = default;
end

%------------------------------------------------------------------
function value = infergetfast(options,name,defaultopt)
%infergetFAST Get infer OPTIONS parameter with no error checking so fast.
%   VAL = infergetFAST(OPTIONS,FIELDNAME,DEFAULTOPTIONS) will get the
%   value of the FIELDNAME from OPTIONS with no error checking or
%   fieldname completion. If the value is [], it gets the value of the
%   FIELDNAME from DEFAULTOPTIONS, another OPTIONS structure which is
%   probably a subset of the options in OPTIONS.
%

if isempty(options)
     value = defaultopt.(name);
     return;
end
% We need to know if name is a valid field of options, but it is faster to use 
% a try-catch than to test if the field exists and if the field name is
% correct. If the options structure is from an older version of the
% toolbox, it could be missing a newer field.
try
    value = options.(name);
catch ME
    value = [];
end

if isempty(value)
    value = defaultopt.(name);
end


