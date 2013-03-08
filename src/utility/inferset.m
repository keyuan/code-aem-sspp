function options = inferset(varargin)
%INFERSET Create/alter inference OPTIONS structure.
%   OPTIONS = inferset('PARAM1',VALUE1,'PARAM2',VALUE2,...) creates an
%   optimization options structure OPTIONS in which the named parameters have
%   the specified values.  Any unspecified parameters are set to [] (parameters
%   with value [] indicate to use the default value for that parameter when
%   OPTIONS is passed to the optimization function). It is sufficient to type
%   only the leading characters that uniquely identify the parameter.  Case is
%   ignored for parameter names.
%   NOTE: For values that are strings, the complete string is required.
%
%   OPTIONS = inferset(OLDOPTS,'PARAM1',VALUE1,...) creates a copy of OLDOPTS
%   with the named parameters altered with the specified values.
%
%   OPTIONS = inferset(OLDOPTS,NEWOPTS) combines an existing options structure
%   OLDOPTS with a new options structure NEWOPTS.  Any parameters in NEWOPTS
%   with non-empty values overwrite the corresponding old parameters in
%   OLDOPTS.
%
%   inferset with no input arguments and no output arguments displays all
%   parameter names and their possible values, with defaults shown in {}
%   when the default is the same for all functions that use that parameter. 
%   Use inferset(OPTIMFUNCTION) to see parameters for a specific function.
%
%   OPTIONS = inferset (with no input arguments) creates an options structure
%   OPTIONS where all the fields are set to [].
%
%   OPTIONS = inferset(OPTIMFUNCTION) creates an options structure with all
%   the parameter names and default values relevant to the optimization
%   function named in OPTIMFUNCTION. For example,
%           inferset('fminbnd')
%   or
%           inferset(@fminbnd)
%   returns an options structure containing all the parameter names and
%   default values relevant to the function 'fminbnd'.
%
%inferset PARAMETERS for inferpp
%Display - Level of display [ off | iter | notify | final ]
%MaxFunEvals - Maximum number of function evaluations allowed
%                     [ positive integer ]
%maxiter - Maximum number of iterations allowed [ positive scalar ]
%tolfun - Termination tolerance on the function value [ positive scalar ]
%tolarg - Termination tolerance on X [ positive scalar ]
%FunValCheck - Check for invalid values, such as NaN or complex, from 
%              user-supplied functions [ {off} | on ]
%OutputFcn - Name(s) of output function [ {[]} | function ] 
%          All output functions are called by the solver after each
%          iteration.
%PlotFcns - Name(s) of plot function [ {[]} | function ]
%          Function(s) used to plot various quantities in every iteration
%
% Note to Optimization Toolbox users:
% To see the parameters for a specific function, check the documentation page 
% for that function. For instance, enter
%   doc fmincon
% to open the reference page for fmincon.
%
% You can also see the options in the Optimization Tool. Enter
%   optimtool
%          
%   Examples:
%     To create an options structure with the default parameters for FZERO
%       options = inferset('fzero');
%     To create an options structure with tolfun equal to 1e-3
%       options = inferset('tolfun',1e-3);
%     To change the Display value of options to 'iter'
%       options = inferset(options,'Display','iter');
%
%   See also OPTIMGET, FZERO, FMINBND, FMINSEARCH, LSQNONNEG.

%   Optimization Toolbox only parameters passed to inferset when the
%   Optimization Toolbox is not on the path now cause a warning (and in a
%   future release an error). To test if the toolbox is on your path, use:
%      ver('optim')

%   Copyright 1984-2009 The MathWorks, Inc.
%   $Revision: 1.34.4.28 $  $Date: 2009/05/18 20:48:13 $


% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    fprintf('                  model: [ ppglm | sspp ]\n');
    fprintf('                  likhd: [ poisson | bernoulli ]\n')
    fprintf('                    cif: [ {} | function]\n');
    fprintf(['                 method: [ em | offvb | onvb | ' ...
            'gibbs | hmc | ...\n']); 
    fprintf(['                           rmhmc | pmmh | lmmh | rbg | '... 
            'irls | cg ]\n'])
    fprintf('                maxiter: [ positive integer ]\n');
    fprintf(['                display: [ off | iter | ', ...
            'notify | final ]\n']);
    fprintf('                 tolfun: [ positive scalar ]\n');
    fprintf('                 tolarg: [ positive scalar ]\n');
    fprintf('                 stadim: [ positive integer ]\n');
    fprintf('               fixparam: [ rho | alpha | sigma | mu | beta ]\n');
    fprintf(['                 intype: [ spike | harmonic | trangle |'... 
            ' box | expdecay | bell ]\n']);
    fprintf('                  estep: [ approxsmoother | varsmoother ]\n');
    fprintf('                 fltopt: [ newton | fixpt | simple ]\n');
    return;
end

% Create a cell array of all the field names
allfields = {'model'; 'likhd'; 'cif'; 'method'; ...
             'maxiter'; 'display'; 'tolfun';'tolarg';'stadim';...
             'fixparam';'intype';'estep';'fltopt'};


% Create a struct of all the fields with all values set to []
% create cell array
structinput = cell(2,length(allfields));
% fields go in first row
structinput(1,:) = allfields';
% []'s go in second row
structinput(2,:) = {[]};
% turn it into correctly ordered comma separated list and call struct
options = struct(structinput{:});

numberargs = nargin; % we might change this value, so assign it
% If we pass in a function name then return the defaults.
if (numberargs==1) && (ischar(varargin{1}) || isa(varargin{1},...
        'function_handle') )
    if ischar(varargin{1})
        funcname = lower(varargin{1});
        if ~exist(funcname,'file')
            error('inferpp:inferset:FcnNotFoundOnPath', ...
                ['No default options available: the function '...
                '%s'' does not exist on the path.'],funcname);
        end
    elseif isa(varargin{1},'function_handle')
        funcname = func2str(varargin{1});
    end
    try 
        optionsfcn = feval(varargin{1},'defaults');
    catch ME
        error('inferpp:inferset:NoDefaultsForFcn', ...
            'No default options available for the function ''%s''.',funcname);
    end
    % The defaults from the optim functions don't include all the fields
    % typically, so run the rest of inferset as if called with
    % inferset(options,optionsfcn)
    % to get all the fields.
    varargin{1} = options;
    varargin{2} = optionsfcn;
    numberargs = 2;
end

Names = allfields;
m = size(Names,1);
names = lower(Names);

i = 1;
while i <= numberargs
    arg = varargin{i};
    if ischar(arg)                         % arg is an option name
        break;
    end
    if ~isempty(arg)                       % [] is a valid options argument
        if ~isa(arg,'struct')
            error('inferpp:inferset:NoParamNameOrStruct',...
                ['Expected argument %d to be a string parameter name ' ...
                'or an options structure\ncreated with inferset.'], i);
        end
        for j = 1:m
            if any(strcmp(fieldnames(arg),Names{j,:}))
                val = arg.(Names{j,:});
            else
                val = [];
            end
            if ~isempty(val)
                if ischar(val)
                    val = lower(deblank(val));
                end
                checkfield(Names{j,:},val);
                options.(Names{j,:}) = val;
            end
        end
    end
    i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(numberargs-i+1,2) ~= 0
    error('inferpp:inferset:ArgNameValueMismatch',...
        'Arguments must occur in name-value pairs.');
end
expectval = 0;                          % start expecting a name, not a value
while i <= numberargs
    arg = varargin{i};

    if ~expectval
        if ~ischar(arg)
            error('inferpp:inferset:InvalidParamName',...
                'Expected argument %d to be a string parameter name.', i);
        end

        lowArg = lower(arg);
        j = strmatch(lowArg,names);
        if isempty(j)                       % if no matches
             error('inferpp:inferset:InvalidParamName',...
                    'Unrecognized parameter name ''%s''.', arg);
        elseif length(j) > 1                % if more than one match
            % Check for any exact matches (in case any names are subsets of others)
            k = strmatch(lowArg,names,'exact');
            if length(k) == 1
                j = k;
            else
                msg = sprintf('Ambiguous parameter name ''%s'' ', arg);
                msg = [msg '(' Names{j(1),:}];
                for k = j(2:length(j))'
                    msg = [msg ', ' Names{k,:}];
                end
                msg = [msg,'.'];
                error('inferpp:inferset:AmbiguousParamName', msg);
            end
        end
        expectval = 1;                      % we expect a value next

    else
        if ischar(arg)
            arg = lower(deblank(arg));
        end
        checkfield(Names{j,:},arg);
        options.(Names{j,:}) = arg;
        expectval = 0;
    end
    i = i + 1;
end

if expectval
    error('inferpp:inferset:NoValueForParam',...
        'Expected value for parameter ''%s''.', arg);
end


%-------------------------------------------------
function checkfield(field,value)
%CHECKFIELD Check validity of structure field contents.
%   CHECKFIELD('field',V,OPTIMTBX) checks the contents of the specified
%   value V to be valid for the field 'field'. 

% empty matrix is always valid
if isempty(value)
    return
end

% See if it is one of the valid inferpp fields.  It may be both an Optim
% and inferpp field, e.g. MaxFunEvals, in which case the inferpp valid
% test may fail and the Optim one may pass.
validfield = 1;
switch field
    case {'tolfun'} % real scalar
        [validvalue, errmsg, errid] = nonNegReal(field,value);
    case {'tolarg'} % real scalar
        % this string is for LSQNONNEG
        [validvalue, errmsg, errid] = nonNegReal(field,value,'10*eps*norm(c,1)*length(c)');
    case {'display'} % several character strings
        [validvalue, errmsg, errid] = displayType(field,value);
    case {'model'} % several character strings
        [validvalue, errmsg, errid] = modeltype(field,value);
    case {'likhd'} % several character strings
        [validvalue, errmsg, errid] = likhdtype(field,value);
    case {'method'} % several character strings
        [validvalue, errmsg, errid] = methodtype(field,value);
    case {'maxiter'} % integer including inf or default string
        % this string is for FMINSEARCH
        [validvalue, errmsg, errid] = nonNegInteger(field,value,'200*numberofvariables');
    case {'cif'}% function
        [validvalue, errmsg, errid] = functionOrCellArray(field,value);
    case {'stadim'} % integer including inf or default string
        % this string is for FMINSEARCH
        [validvalue, errmsg, errid] = nonNegInteger(field,value,'200*numberofvariables');
    case {'fixparam'} % several character strings
        [validvalue, errmsg, errid] = fixparamtype(field,value);
    case {'intype'} % several character strings       
        [validvalue, errmsg, errid] = intypetype(field,value);
    case {'estep'} % several character strings       
        [validvalue, errmsg, errid] = esteptype(field,value);
    case {'fltopt'} % several character strings       
        [validvalue, errmsg, errid] = fltopttype(field,value);
    otherwise
        validfield = 0;  
        validvalue = 0;
        errmsg = sprintf('Unrecognized parameter name ''%s''.', field);
        errid = 'inferpp:inferset:checkfield:InvalidParamName';
end

if validvalue 
    return;
elseif  validfield  
    % Throw the inferpp invalid value error
    error(errid, errmsg);
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = nonNegReal(field,value,string)
% Any nonnegative real scalar or sometimes a special string
valid =  isreal(value) && isscalar(value) && (value >= 0) ;
if nargin > 2
    valid = valid || isequal(value,string);
end

if ~valid
    if ischar(value)
        errid = 'inferpp:funfun:inferset:NonNegReal:negativeNum';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: ... must be a real non-negative scalar (not a string).',field);
    else
        errid = 'inferpp:funfun:inferset:NonNegReal:negativeNum';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-negative scalar.',field);
    end
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = nonNegInteger(field,value,string)
% Any nonnegative real integer scalar or sometimes a special string
valid =  isreal(value) && isscalar(value) && (value >= 0) && value == floor(value) ;
if nargin > 2
    valid = valid || isequal(value,string);
end
if ~valid
    if ischar(value)
        errid = 'inferpp:funfun:inferset:nonNegInteger:notANonNegInteger';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-negative integer (not a string).',field);
    else
        errid = 'inferpp:funfun:inferset:nonNegInteger:notANonNegInteger';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-negative integer.',field);
    end
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = displayType(field,value)
% One of these strings: on, off, none, iter, final, notify
valid =  ischar(value) && any(strcmp(value, ...
    {'notify';'off';'iter';'final'}));
if ~valid
    errid = 'inferpp:funfun:inferset:displayType:notADisplayType';
    errmsg = sprintf(['Invalid value for OPTIONS parameter %s: must be ''off'',''on'',''iter'',\n', ...
        '''iter-detailed'',''notify'',''notify-detailed'',''final'', or ''final-detailed''.'],field);
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = onOffType(field,value)
% One of these strings: on, off
valid =  ischar(value) && any(strcmp(value,{'on';'off'}));
if ~valid
    errid = 'inferpp:funfun:inferset:onOffType:notOnOffType';
    errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''off'' or ''on''.',field);
else
    errid = '';
    errmsg = '';
end

%--------------------------------------------------------------------------------

function [valid, errmsg, errid] = modeltype(field,value)
% One of these strings: on, off
valid =  ischar(value) && any(strcmp(value,{'ppglm';'sspp'}));
if ~valid
    errid = 'inferpp:funfun:inferset:onOffType:notOnOffType';
    errmsg = sprintf(['Invalid value for OPTIONS parameter %s:' ...
                     'must be ''sspp'' or ''ppglm''.'],field);
else
    errid = '';
    errmsg = '';
end

%--------------------------------------------------------------------------

function [valid, errmsg, errid] = likhdtype(field,value)
% One of these strings: on, off
valid =  ischar(value) && any(strcmp(value,{'poisson';'bernoulli'}));
if ~valid
    errid = 'inferpp:funfun:inferset:onOffType:notOnOffType';
    errmsg = sprintf(['Invalid value for OPTIONS parameter %s:' ...
                     'must be ''poisson'' or ''bernoulli''.'],field);
else
    errid = '';
    errmsg = '';
end

%--------------------------------------------------------------------------

function [valid, errmsg, errid] = methodtype(field,value)
% One of these strings: on, off
valid =  ischar(value) && any(strcmp(value,{'em';'offvb';'onvb';...
    'gibbs';'hmc';'rmhmc';'pmmh';'lmmh';'rbg';'irls';'cg'}));
if ~valid
    errid = 'inferpp:funfun:inferset:onOffType:notOnOffType';
    errmsg = sprintf(['Invalid value for OPTIONS parameter %s:' ...
                     'must be ''em'', ''offvb'', ''onvb''' ...
    ', ''gibbs'', ''hmc'', ''rmhmc'', ''pmmh'', ''lmmh'', ''rbg'', ''irls'''...
    ' or ''cg''.'],field);
else
    errid = '';
    errmsg = '';
end

%--------------------------------------------------------------------------

function [valid, errmsg, errid] = fixparamtype(field,value)
% One of these strings: on, off
for ii = 1:size(value,1)
valid = ischar(value{ii}) && any(strcmp(value{ii},{'rho';'alpha';'sigma';'mu';'beta'}));
if ~valid
    errid = 'inferpp:funfun:inferset:onOffType:notOnOffType';
    errmsg = sprintf(['Invalid value for OPTIONS parameter %s:' ...
                     'must be ''rho'', ''alpha'', ''sigma'', ' ...
                     '''mu'' or ''beta''.'],field);
else
    errid = '';
    errmsg = '';
end
end
%--------------------------------------------------------------------------

function [valid, errmsg, errid] = intypetype(field,value)
% One of these strings: on, off
valid =  ischar(value) && any(strcmp(value,{'spike';'harmonic';'trangle';...
    'box';'expdecay';'bell';'rbf';'sigmoid';'none'}));
if ~valid
    errid = 'inferpp:funfun:inferset:onOffType:notOnOffType';
    errmsg = sprintf(['Invalid value for OPTIONS parameter %s:' ...
                     'must be ''spike'', ''harmonic''' ...
    ', ''trangle'', ''box'', ''expdecay'', ''bell'', ''rbf'', ''sigmoid'', '...
    ' ''none''. '],field);
else
    errid = '';
    errmsg = '';
end

%----------------------------------------------------------------

function [valid, errmsg, errid] = esteptype(field,value)
% One of these strings: on, off
valid =  ischar(value) && any(strcmp(value,{'approxsmoother';'varsmoother'}));
if ~valid
    errid = 'inferpp:funfun:inferset:onOffType:notOnOffType';
    errmsg = sprintf(['Invalid value for OPTIONS parameter %s:' ...
                     'must be ''approxsmoother'' or ''varsmoother''.'],field);
else
    errid = '';
    errmsg = '';
end

%--------------------------------------------------------------------------

function [valid, errmsg, errid] = fltopttype(field,value)
% One of these strings: on, off
valid =  ischar(value) && any(strcmp(value,{'newton';'fixpt';'simple'}));
if ~valid
    errid = 'inferpp:funfun:inferset:onOffType:notOnOffType';
    errmsg = sprintf(['Invalid value for OPTIONS parameter %s:' ...
                     'must be ''newton'', ''fixpt'' or ''simple''.'],field);
else
    errid = '';
    errmsg = '';
end

%--------------------------------------------------------------------------


function [valid, errmsg, errid] = functionOrCellArray(field,value)
% Any function handle, string or cell array of functions 
valid =  ischar(value) || isa(value, 'function_handle') || iscell(value);
if ~valid
    errid = 'inferpp:funfun:inferset:functionOrCellArray:notAFunctionOrCellArray';
    errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a function or a cell array of functions.',field);
else
    errid = '';
    errmsg = '';
end
%--------------------------------------------------------------------------------


