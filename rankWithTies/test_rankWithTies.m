% test_rankWithTies - Test suite
%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   This code provides the test suite.  It can be run through the MOxUnit
%   unit testing framework.
%
% PURPOSE:
%   This function provides the unit tests for cconv.m  New bug reports for the
%   code should not be closed with ensuring that a test in the suite both fails
%   before the fix, and passes after it.
%
%-------------------------------------------------------------------------------
%}
function test_rankWithTies()

try % assignment of 'localfunctions' is necessary in Matlab >= 2016
    test_functions = localfunctions();
catch % no problem; early Matlab versions can use initTestSuite fine
end
initTestSuite;

end



% Test that all of the described syntax is supported
function test_syntaxSupport()

ACell = {{[pi 1 1 2 3 3 ]}};
methodCell = {{'min'}, {'max'}, {'random'}, {'stable'}, {@(rMax,nMatch)rMax} };
sortOptsCell = {{'descend'},{}};

 % Iterate ensuring every permutation of calling syntaxes execute without error
 argSize = [numel(ACell) numel(methodCell) numel(sortOptsCell)];
 
 for i = 1:prod(argSize)
    argInd = ind2subVec(argSize, i, true);
    args = [...
        ACell{argInd(1)},       ...
        methodCell{argInd(2)},  ...
        sortOptsCell{argInd(3)} ];
    try
        % This shouldn't throw
        [R,B,I] = rankWithTies( args{:} );
        
    catch err
        % Tell users which case failed
        fprintf('Unsupported argument set: %d\n', i);
        
        % Give the error
        rethrow(err);
    end
 end


end


% Test type support
function test_typeSupport()

A = 'fooBar...';
ACell = {{double(A)}, {single(A)}, {uint8(A)}, {char(A)}, ...
    {logical([1 1 0 1 0 0])}};
methodCell = {{'min'}, {'max'}, {'random'}, {'stable'}, {@(rMax,nMatch)rMax} };
sortOptsCell = {{'descend'},{}};

 % Iterate ensuring every permutation of calling syntaxes execute without error
 argSize = [numel(ACell) numel(methodCell) numel(sortOptsCell)];
 
 for i = 1:prod(argSize)
    argInd = ind2subVec(argSize, i, true);
    args = [...
        ACell{argInd(1)},       ...
        methodCell{argInd(2)},  ...
        sortOptsCell{argInd(3)} ];
    
    try
        % This shouldn't throw
        [R,B,I] = rankWithTies( args{:} );
        
        % Verify the return types
        assertEqual( class(R), 'double' );
        assertEqual( class(B), class(ACell{argInd(1)}{:}) );
        assertEqual( class(I), 'double' );
        
    catch err
        % Tell users which case failed
        fprintf('Unsupported argument set: %d\n', i);
        
        % Give the error
        rethrow(err);
    end
 end

end


% Test dimension support
function test_dimensionSupport

A = 'fooBar...';
nA = numel(A);

for dimInd = 1:4
    % Setup A
    aSize = ones(1,max(2,dimInd));
    aSize(dimInd) = nA;
    A = reshape(A,aSize);
    
    % Confirm the dimensions match
    [R,B,I] = rankWithTies( A, 'min');
    assertEqual( size(R), aSize );
    assertEqual( size(B), aSize );
    assertEqual( size(I), aSize );
    
end

end


% Test permutation properties
function test_permutationProperties

A = 'fooBar...';
method = {'min','max','random','stable',@(rMax,nMatch) (rMax-nMatch):rMax};

for i = 1:numel(method)
    [R,B,I] = rankWithTies(A, method{i});
    assertEqual( B(R), A);
    assertEqual( A(I), B);
end

end


% Test method function equivalents
% Note: We need to carefully manage the random number generator seed
function test_methodFunctionalEquivalents

method = {'min','max','random','stable'};
equivalent = {...
    @(rMax,nMatch) rMax-nMatch, ...
    @(rMax,nMatch) rMax, ...
    @(rMax,nMatch) randperm(nMatch+1) + (rMax-nMatch-1), ...
    @(rMax,nMatch) (rMax-nMatch):rMax };
A = 'fooBar...';

for i = 1:numel(method)
    % Note the random number generator state
    rngState = rng();
    
    % Rank using method
    [R1,B1,I1] = rankWithTies( A, method{i});
    
    % Reset the random number generator state
    rng(rngState);
    
    % Rank using equivalent
    [R2,B2,I2] = rankWithTies( A, equivalent{i});
    
    % Compare
    assertEqual( R1, R2);
    assertEqual( B1, B2);
    assertEqual( I1, I2);
end

end


% Test sort options
% Note: We conspire to ensure vectorization is not required via dim
function test_sortOptions

% Default Values
A = 'fooBar...';
nA = numel(A);
method = 'min';
dimOpts = {{1}, {2}, {}};
modeOpts = {{'ascend'}, {'descend'}, {}};

% Get some true answers associated with each sort mode
expectedB = {sort(A,'ascend'), sort(A,'descend'), sort(A)};

% Iterate ensuring every permutation of calling syntaxes execute without error
 argSize = [numel(dimOpts) numel(modeOpts)];
 
 for i = 1:prod(argSize)
    argInd = ind2subVec(argSize, i, true);
    sortArgs = [ dimOpts{argInd(1)}, modeOpts{argInd(2)} ];
    
    % Be sure A is in the dimension of dimOpts (no vectorization)
    if cellfun(@isempty, dimOpts(argInd(1)))
        d = 1;
    else
        d = dimOpts{argInd(1)}{1};
    end
    ASize = ones(1,max(2,d));
    ASize(d) = nA;
    A = reshape(A,ASize);
    
    % Verify the call completes
    [R,B,I] = rankWithTies(A, method, sortArgs{:});
    
    % Verify permutation properties
    assertEqual( B(R), A);
    assertEqual( A(I), B);
    
    % Verify the expected sort result
    assertEqual( B(:), expectedB{argInd(2)}(:) );
    
 end

end


% Test vectorization support
function test_vectorization

% Default Values
method = 'min';


% Vectorization along dimension 1
ASize = [5 6];
A = randn(ASize);
for i = ASize(2):-1:1
    [RTrue(:,i), BTrue(:,i), ITrue(:,i)] = rankWithTies(A(:,i), method);
end

[R,B,I] = rankWithTies(A, method);      % Implicitly on dimension 1
assertEqual(R, RTrue);
assertEqual(B, BTrue);
assertEqual(I, ITrue);

[R,B,I] = rankWithTies(A, method, 1);   % Explicitly on dimension 1
assertEqual(R, RTrue);
assertEqual(B, BTrue);
assertEqual(I, ITrue);


% Vectorization along dimension 2
for i = ASize(1):-1:1
    [RTrue(i,:), BTrue(i,:), ITrue(i,:)] = rankWithTies(A(i,:), method);
end

[R,B,I] = rankWithTies(A, method, 2);   % Explicitly on dimension 2
assertEqual(R, RTrue);
assertEqual(B, BTrue);
assertEqual(I, ITrue);

% Trivial vectorization
[R,B,I] = rankWithTies(A, method, 3);
assertEqual(R, ones(ASize));
assertEqual(B, A);
assertEqual(I, ones(ASize));

end


% Test edge cases
function test_edgeCases

method = 'max';
A = [];
ASize = size(A);

[R,B,I] = rankWithTies(A, method);
assertEqual( size(R), ASize);
assertEqual( size(B), ASize);
assertEqual( size(I), ASize);

[R,B,I] = rankWithTies(A, method, 1);
assertEqual( size(R), ASize);
assertEqual( size(B), ASize);
assertEqual( size(I), ASize);

[R,B,I] = rankWithTies(A, method, 2);
assertEqual( size(R), ASize);
assertEqual( size(B), ASize);
assertEqual( size(I), ASize);

[R,B,I] = rankWithTies(A, method, 3);
assertEqual( size(R), ASize);
assertEqual( size(B), ASize);
assertEqual( size(I), ASize);

A = 1;
ASize = size(A);

[R,B,I] = rankWithTies(A, method);
assertEqual(R, 1);
assertEqual(B, 1);
assertEqual(I, 1);

[R,B,I] = rankWithTies(A, method, 1);
assertEqual(R, 1);
assertEqual(B, 1);
assertEqual(I, 1);

[R,B,I] = rankWithTies(A, method, 2);
assertEqual(R, 1);
assertEqual(B, 1);
assertEqual(I, 1);

[R,B,I] = rankWithTies(A, method, 3);
assertEqual(R, 1);
assertEqual(B, 1);
assertEqual(I, 1);


end

