function z = plus( x, y )
    %PLUS Addition of two TT/MPS tensors.
    %   Z = PLUS(X,Y) adds two TT/MPS tensors. The rank of the resulting
    %   tensor is 2*R.
    %
    %   See also MINUS, UMINUS.
    
    %   TTeMPS Toolbox. 
    %   Michael Steinlechner, 2013-2016
    %   Questions and contact: michael.steinlechner@epfl.ch
    %   BSD 2-clause license, see LICENSE.txt
    
    % add sanity check...
    rx = x.rank;
    ry = y.rank;
    nx = x.size;

    z = TTeMPS( cell(1, x.order) );
        
    % first core:
    p = size(x.U{1},4);
    tmp = zeros( 1, nx(1), rx(2)+ry(2), p );
    tmp( 1, :, 1:rx(2), : ) = x.U{1};
    tmp( 1, :, rx(2)+1:end, : ) = y.U{1};
    z.U{1} = tmp;

    % central cores:
    for i = 2:x.order-1
        % possibility of block format:
        p = size(x.U{i},4);
        tmp = zeros( rx(i)+ry(i), nx(i), rx(i+1)+ry(i+1), p);
        tmp( 1:rx(i), :, 1:rx(i+1), :) = x.U{i};
        tmp( rx(i)+1:end, :, rx(i+1)+1:end, :) = y.U{i};
        z.U{i} = tmp;
    end

    % last core:
    p = size(x.U{end},4);
    tmp = zeros( rx(end-1)+ry(end-1), nx(end), 1, p );
    tmp( 1:rx(end-1), :, 1, : ) = x.U{end};
    tmp( rx(end-1)+1:end, :, 1, : ) = y.U{end};
    z.U{end} = tmp;
end
