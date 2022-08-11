%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   z = monotonicRegression(x, target, floor)
%
% PURPOSE:
%   If floor is true (default) floors x to target if x > target, else leaves
%       it unchanged. If floor is false rounds x up to target if x < target
%
% INPUT:
%   x        - Matrix of numeric values
%   target   - Value to floor or ceil to. Same dimension of x or single value
%   floor    - Boolean toggling floor or ceil modes
%
% OUTPUT:
%   z        - min(0, (x - target)) + target or max(0, (x - target)) + target;
%-------------------------------------------------------------------------------
%}
function z = roundTo(x, target, varargin)

floor = parseInputs(varargin{:});

if floor
    z = min(0, (x - target)) + target;
else
    z = max(0, (x - target)) + target;
end


end

function floor = parseInputs(varargin)

    switch nargin 
        case 0
           floor = true;
        case 1
           floor = varargin{1};
        otherwise
            error([mfilename,':invalidInput'],...
                'Unexpected number of inputs');
    end
end
