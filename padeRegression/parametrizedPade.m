%{
%-------------------------------------------------------------------------------
% SYNTAX:
%   [S,<G>] = padeReformed(<x>,cm)
%
% PURPOSE:
%   Computes the Pade approximant of order [2/1] defined by its starting point
%       (x0, y0), midpoint (xc, yc), and endpoint (xf, yf)
%
% INPUTS:
%   x  - value of x to evaluate the spline at 
%   xc - x value of the control point
%   yc - y value of the control point
%   x0 - x value of the starting point of the spline
%   y0 - y value of the starting point of the spline
%   xf - x value of the ending point of the spline
%   yf - y value of the ending point of the spline
%
% OUTPUTS:
%   S       - A vector containing the values of the splines evaluated at x if x
%             is in that splines respective domain.
%   G       - A vector containing the values of the gradients of the splines
%             evaluated at x if x is in that splines respective domain.
%
% NOTE:
%   The spline calculated is constrained to have a continuous second derivative
%       and thus be one Pade approximant 
%
%-------------------------------------------------------------------------------
%}
function [y, yprime] = parametrizedPade(x, xc, yc, varargin)

[x0, y0, xf, yf] = parseInputs(varargin{:});

y = (x.*x0.*y0.*yc - x.*x0.*y0.*yf - x.*xc.*y0.*yc + x0.*xc.*y0.*yf - x0.*xf.*y0.*yc + ...
x.*xf.*y0.*yf + x.*xc.*yc.*yf - x0.*xc.*yc.*yf + xc.*xf.*y0.*yc - x.*xf.*yc.*yf+ ...
x0.*xf.*yc.*yf - xc.*xf.*y0.*yf) ./ (x.*x0.*yc - x.*xc.*y0 + x0.*xc.*y0 - x.*x0.*yf + ...
x.*xf.*y0 - x0.*xf.*y0 - x0.*xc.*yc + x.*xc.*yf - x.*xf.*yc + x0.*xf.*yf + xc.*xf.*yc ...
- xc.*xf.*yf);

% cross checked y,yprime against oneNodePadeSplines and see agreement
yprime = ((xc - x0) * (xf - x0) * (xf - xc) * (yc - y0) * (yf - y0) * (yf - yc)) ...
            ./ ((((xc - x0)*yf - (xf - x0)*yc + (xf - xc)*y0)*x ...
            - (xc - x0)*xf*yf + xc*yc*(xf - x0) - x0*y0*(xf - xc)).^2);

end

function [x0, y0, xf, yf] = parseInputs(varargin)
    
    switch numel(varargin)
        case 0 
            x0 = 0;
            y0 = 0;
            xf = 1;
            yf = 1;
        case 4
            x0 = varargin{1};
            y0 = varargin{2};
            xf = varargin{3};
            yf = varargin{4};
        otherwise 
            error([mfilename,':invalidInput'],...
            'Unexpected number of inputs');
    end

end