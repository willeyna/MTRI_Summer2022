% rankWithTies - Ranking with tie-break control
%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   [R,B,I] = rankWithTies(A,method,<sortOpts>)
%
% PURPOSE:
%   This function performs ranking with additional control of tie-breaks.  
%   Ranking is fundamentally a sorting operator where one reports the relative
%   position of each element within the sorted list.  When ties are present,
%   these relative positions become ambiguous, and this function controls how 
%   those positions get reported.
% 
%   The outputs R and I are indices that describe permutation operators with the
%   properties A = B(R) and B = A(I) for vector-valued inputs.  R and I will not
%   represent inverse permutation operators when ties are present if one does
%   not ensure unique(R) = unique(I).  For example, the 'min' and 'max' methods
%   below will lead to R and I representing permutations that are not inverses
%   of one another when ties are present.
%  
% INPUT:
%   A           - Array to be sorted
%
%   method      - Method for breaking ties.  Common defaults are provided and 
%                 one can also define an arbitrary function that acts on the
%                 stable ranks
%       'min'       - Smallest index               
%       'max'       - Largest index
%       'random'    - Random permutation index
%       'stable'    - Order of occurance in A
%
%                 A custom function defined on the stable ranks given by
%                 (rMax-nMatch):rMax.  This function must have to prototype
%                   ROut(ind) = methodFun(rMax,nMatch)
%                 where rMax is the maximum stable rank and nMatch is the
%                 number of matches within the set B( (rMax-nMatch):rMax ).
%                 For the aforementioned defaults this function is:
%       'min'       - @(rMax,nMatch) rMax-nMatch
%       'max'       - @(rMax,nMatch) rMax
%       'random'    - @(rMax,nMatch) randperm(nMatch+1) + (rMax-nMatch-1)
%       'stable'    - @(rMax,nMatch) (rMax-nMatch):rMax 
%                 
%
%   sortOpts    - Options to built-in sort call
% 
% OUTPUT:
%   R       - Ranks.  Reverse sort index s.t. A = B(R)
%   B       - Sort result
%   I       - Sort index with ties broken via method s.t. B = A(I)
%
% ASSUMPTIONS: 
%   All input variables are of the correct type, valid (if applicable),
%   and given in the correct order.
%
% NOTES:
%   Vectorization doesn't explicitly parallelize replaceDuplicateRanks 
%
% SEE ALSO:
%   datafun/sort
%-------------------------------------------------------------------------------
%}
function [R,B,I] = rankWithTies(A,method,varargin)

% Interpret cell array sortOpts as a cell of options.  This nicety works because
% sort does not accept options of type cell.
if nargin==3 && iscell(varargin{1})
    varargin = varargin{1};
end

% Determine if vectorization is being used
vecDim = getVectorizationDim(A,varargin{:});

% Sort
[B,I] = sort(A,varargin{:});

% Get the stable ranks (reverse indices)
R = getStableRanks(I,vecDim);

% Get the duplicate rank function
switch class(method)
    case 'char'
        % Provide the defaults
        switch lower(method)
            case 'min'
                % min( R(I( (rMax-nMatch):rMax )) )
                fun = @(rMax,nMatch) rMax-nMatch;
                
            case 'max'
                % max( R(I( (rMax-nMatch):rMax )) )
                fun = @(rMax,nMatch) rMax;
                
            case {'random','rand'}
                % Random permutation of R(I( (rMax-nMatch):rMax  ))
                fun = @(rMax,nMatch) randperm(nMatch+1) + (rMax-nMatch-1);
                
            case 'stable'
                % Nothing left to do
                return
                
            otherwise
                error([mfilename(),':unknownMethod'],...
                    'Unknown tie-break method: %s', ...
                    strtrim(evalc('disp(method)')) );
        end
        
    case 'function_handle'
        % User provided function
        fun = method;
        
    otherwise
        error([mfilename(),':invalidMethodType'],...
                    'Invalid method type: %s', class(method) );
end

% Replace duplicate stable ranks according to fun
if vecDim
    R = replaceDuplicateRanksVec(B, I, R, fun, vecDim);
else
    R = replaceDuplicateRanks(B, I, R, fun );
end


end


%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   vecDim = getVectorizationDim(A,varargin)
%
% PURPOSE:
%   Returns the working dimension if vectorization is being used.
%  
% INPUT:
%   See main function
% 
% OUTPUT:
%   vecDim      - Vectorization dimension or [] if no vectorization
%
% ASSUMPTIONS: 
%   All input variables are of the correct type, valid (if applicable),
%   and given in the correct order.
%-------------------------------------------------------------------------------
%}
function vecDim = getVectorizationDim(A,varargin)

% Default Values
ASize = size(A);
    
if numel(varargin) && isnumeric(varargin{1})
    % User provided an explicit dimension to sort
    vecDim = varargin{1};
    
    % Check for vectorization
    if vecDim <= numel(ASize)
        ASize(vecDim) = [];
    end
    if prod(ASize)==1
        vecDim = [];
    end
    
elseif sum(ASize>1)>1
    vecDim = 1;
    
else
    vecDim = [];
end

end


%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   R = getStableRanks(I,vecDim)
%
% PURPOSE:
%   Returns the stable ranks associated with I.  This is the reverse permutation
%   indices s.t. R(I) are ranks of elements that are sorted by I.
%  
% INPUT:
%   I           - Sort index s.t. B = A(I) sorts A
%   vecDim      - Vectorization dimension or [] if no vectorization
% 
% OUTPUT:
%   R           - Stable ranks s.t. A = B(R)
%   IFull       - I described in terms of the full vectorized problem
%
% ASSUMPTIONS: 
%   All input variables are of the correct type, valid (if applicable),
%   and given in the correct order.
%-------------------------------------------------------------------------------
%}
function [R,I] = getStableRanks(I,vecDim)

% Default Values
ISize = size(I);

if isempty(vecDim)
    % Get stable reverse indices
    N = numel(I);
    R(I) = 1:N;
    
    % Ensure R is of the correct size
    if ISize(2)~=N || ISize(2)==0
        R = reshape(R,ISize);
    end
    
else    
    % Trivial vectorization
    if vecDim > numel(ISize)
        R = I;
        return
    end
    
    % Get the vectorized assignments "1:N"
    nDim = numel(ISize);
    for i = 1:nDim
        if i==vecDim
            ind{i} = ':';
        else
            ind{i} = ones(ISize(i),1);
        end
    end
    RSize = ones(1,nDim);
    RSize(vecDim) = ISize(vecDim);
    R = reshape(1:ISize(vecDim), RSize);
    R = R(ind{:});
    
    % Compute the delta along vecDim
    delta = prod(ISize(1:vecDim-1));
    
    % Convert I to the linear index associated with the full problem
    I = reshape(1:prod(ISize),ISize) + (I-R) * delta;
    
    % Transform R into the desired reverse mapping
    R(I) = R;
    
end

end


%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   R = replaceDuplicateRanks(B,I,R,fun)
%
% PURPOSE:
%   Replace duplicate ranks with a function of their stable ordering.
%  
% INPUT:
%   B       - Sorted output
%   I       - Forward sort index s.t. B = A(I)
%   R       - Stable reverse sort index s.t. A = B(R)
%   fun     - Function mapping stable ranks to their desired output where the
%             stable ranks are given by (rMax-nMatch):rMax
%               ROut(ind) = fun(rMax,nMatch)
% 
% OUTPUT:
%   R       - Updated reverse sort index (ranks) s.t. A = B(R)
%
% ASSUMPTIONS: 
%   All input variables are of the correct type, valid (if applicable),
%   and given in the correct order.
%-------------------------------------------------------------------------------
%}
function R = replaceDuplicateRanks(B,I,R,fun)

% Default Values
N = numel(B);
matchCount = 0;

for i = 1:N-1
    if B(i)==B(i+1)
        matchCount = matchCount + 1;
        
    elseif matchCount
        % Resolve the matched indices
        % ind = (i-matchCount):i;
        % R(I(ind)) = fun( R(I(ind)) );
        R(I( (i-matchCount):i )) = fun(i,matchCount);
        
        matchCount = 0;
    end

end

% Resolve remaining matches
if matchCount
    R(I( (N-matchCount):N )) = fun(N,matchCount);
end

end



%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   R = replaceDuplicateRanksVec(B,I,R,fun,vecDim)
%
% PURPOSE:
%   Vectorize replaceDuplicateRanks.  The most efficient path is not yet clear.
%  
% INPUT:
%   B       - Sorted output
%   I       - Forward sort index s.t. B = A(I)
%   R       - Stable reverse sort index s.t. A = B(R)
%   fun     - Function mapping stable ranks to their desired output where the
%             stable ranks are given by (rMax-nMatch):rMax
%               ROut(ind) = fun(rMax,nMatch)
%   vecDim  - Dimension of vectorization
%
% OUTPUT:
%   R       - Updated reverse sort index (ranks) s.t. A = B(R)
%
% ASSUMPTIONS: 
%   All input variables are of the correct type, valid (if applicable),
%   and given in the correct order.
%-------------------------------------------------------------------------------
%}
function R = replaceDuplicateRanksVec(B,I,R,fun,vecDim)

% Default Values
BSize = size(B);
nDim = numel(BSize);

if numel(BSize)<vecDim
    % Trivial vectorization
    R = replaceDuplicateRanks(B,I,R,fun);
    return
elseif numel(B) == 0
    % Trivial vectorization
    R = zeros(BSize);
    return
end

% Place the working dimension up front
if prod(BSize(1:vecDim-1)) > 1
    order = 0:nDim;
    order(vecDim+1) = [];
    order(1) = vecDim;
    
    B = permute(B,order);
    I = permute(I,order);
    R = permute(R,order);
else
    order = [];
end

% Iterate over each slice
N = numel(B)/BSize(vecDim);
for i = 1:N
    R(:,i) = replaceDuplicateRanks(B(:,i),I(:,i),R(:,i),fun);
end

% Inverse permute
if order
    R = ipermute(R,order);
end

end