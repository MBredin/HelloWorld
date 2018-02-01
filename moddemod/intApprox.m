function y = intApprox(x,fs,method)
%==========================================================================
% Call Syntax: y = intApprox(x,fs,method)
%
% Description: This function evaluates an approximate computation of an 
%              integral using various numerical techniques.
% 
% Input Arguments:
%   Name: x
%   Type:  column vector, or matrix of column vectors 
%   Description: signal(s) to be intergrated 
%
%   Name: fs (optional)
%   Type: scalar
%   Description: sampling frequency [default = 1]
%
%   Name: method (optional)
%   Type: string
%   Description: numerical integration method:
%                   'left'   - use the left rectangle rule
%                   'right'  - use the right rectangle rule
%                   'center' - use the midpoint rule
%                   'trapz'  - use the trapazoidal rule [default]  
%                   'simps'  - use Simpson's rule
%
% Output Arguments:
%   Name: y
%   Type: column vector, or matrix of column vectors 
%   Description:integral estimate       
%
%
%--------------------------------------------------------------------------
% If you use these files please cite the following:
%
%       @article{HSA2017,
%           title={The Hilbert Spectrum: A General Framework for Time-Frequency Analysis},
%           author={Sandoval, S. and De~Leon, P.~L.~},
%           journal={{IEEE Trans.~Signal Process.}},
%           year = {\noop{2017}in review},  }
%
%--------------------------------------------------------------------------
%
% References:
%
%
% Notes: See DEMO_intApprox.m
%
%
% Function Dependencies:  
%
%--------------------------------------------------------------------------
% Author: Steven Sandoval
%--------------------------------------------------------------------------
% Creation Date: May 2017
%
% Revision History:  9-11-2017 S. Terrazas - updated function call
%
%==========================================================================

%------------------
% Check valid input
%------------------
if nargin<3
    method = 'trapz';
end

if nargin<2
    fs = 1;
end

if not(strcmp(method,'left')||strcmp(method,'right')||strcmp(method,'center')||strcmp(method,'trapz')||strcmp(method,'simps'))
    error('Error (intApprox): invalid integration method.')
end

%DEAL WITH MATRIX INPUT
rows=size(x,1);
cols=size(x,2);
if cols>1
    y=zeros(rows, cols);
    for k=1:size(x,2)
        [y(:,k)]=intApprox(x(:,k), fs,method);
    end
    return
end

%-----------
% Initialize
%-----------


%-----
% Main
%-----

switch method
    case 'left'
        y = cumsum(x);
    case 'right'
        yi1 = cumsum(x);
        y = [0;yi1(1:end-1)];
    case 'center'
        yi1 = cumsum(x);
        yi2 = [0;yi1(1:end-1)];
        y = 1/2*(yi1+yi2);
    case 'trapz'
        y = cumtrapz(x);
    case 'simps'
        y = cumsimps(x);
end
y = y./fs;

end



%==============================SUB-ROUTINES=================================
% This sub-routine was taken from : http://www.biomecardio.com/matlab/cumsimps.html

function z = cumsimps(x,y,dim)

%CUMSIMPS Cumulative Simpson's numerical integration.
%   Z = CUMSIMPS(Y) computes an approximation of the cumulative integral of
%   Y via the Simpson's method (with unit spacing). To compute the integral
%   for spacing different from one, multiply Z by the spacing increment.
%
%   For vectors, CUMSIMPS(Y) is a vector containing the cumulative integral
%   of Y. For matrices, CUMSIMPS(Y) is a matrix the same size as X with the
%   cumulative integral over each column. For N-D arrays, CUMSIMPS(Y) works
%   along the first non-singleton dimension.
%
%   Z = CUMSIMPS(X,Y) computes the cumulative integral of Y with respect to
%   X using Simpson's integration. X and Y must be vectors of the same
%   length, or X must be a column vector and Y an array whose first
%   non-singleton dimension is length(X). CUMSIMPS operates across this
%   dimension.
%
%   Z = CUMSIMPS(X,Y,DIM) or CUMSIMPS(Y,DIM) integrates along dimension DIM
%   of Y. The length of X must be the same as size(Y,DIM).
%
%   Examples:
%   --------
%   % The integral of cos(x) is sin(x)
%   % Let us compare CUMTRAPZ and CUMSIMPS
%   x = linspace(0,2*pi,20);
%   y = cos(x);
%   yt = cumtrapz(x,y);
%   ys = cumsimps(x,y);
%   % RMSE: root mean squared errors:
%   sqrt(mean((yt-sin(x)).^2)) % returns 0.0063
%   sqrt(mean((ys-sin(x)).^2)) % returns 1.2309e-004
%
%   If Y = [0 1 2; 3 4 5; 6 7 8]
%   then cumsimps(Y,1) is [0   0   0    and cumsimps(Y,2) is [0 0.5 2
%                          1.5 2.5 3.5                        0 3.5 8
%                          6   8   10]                        0 6.5 14]
%
%
%   Class support for inputs X,Y:
%      float: double, single
%
%   -- Damien Garcia -- 09/2007, revised 11/2009
%   directly adapted from CUMTRAPZ
%
%   See also CUMSUM, CUMTRAPZ, SIMPS.

%   Make sure x and y are column vectors, or y is a matrix.
perm = []; nshifts = 0;
if nargin == 3, % cumsimps(x,y,dim)
    perm = [dim:max(length(size(y)),dim) 1:dim-1];
    y = permute(y,perm);
    [m,n] = size(y);
elseif nargin==2 && isequal(size(y),[1 1]) % cumsimps(y,dim)
    dim = y; y = x;
    perm = [dim:max(length(size(y)),dim) 1:dim-1];
    y = permute(y,perm);
    [m,n] = size(y);
    x = 1:m;
else
    if nargin < 2, y = x; end
    [y,nshifts] = shiftdim(y);
    [m,n] = size(y);
    if nargin < 2, x = 1:m; end
end
x = x(:);
if length(x) ~= m
    error('MATLAB:cumsimps:LengthXMismatchY',...
        'length(x) must equal length of first non-singleton dim of y.');
end
if m<3
    error('MATLAB:cumsimps:LengthX',...
        'length(x) must be at least 3.');
end

%
z = zeros(size(y),class(y));

dx = repmat(diff(x,1,1),1,n);
dx1 = dx(1:end-1,:);
dx2 = dx(2:end,:);

alpha = (dx1+dx2)./dx1/6;
a0 = alpha.*(2*dx1-dx2);
a1 = alpha.*(dx1+dx2).^2./dx2;
a2 = alpha.*dx1./dx2.*(2*dx2-dx1);

% First cumulative value
state0 = warning('query','MATLAB:nearlySingularMatrix');
state0 = state0.state;
warning('off','MATLAB:nearlySingularMatrix')
C = vander(x(1:3))\y(1:3,:);
z(2,:) = C(1,:).*(x(2,:).^3-x(1,:).^3)/3 +...
    C(2,:).*(x(2,:).^2-x(1,:).^2)/2 +...
    C(3,:).*dx(1,:);
warning(state0,'MATLAB:nearlySingularMatrix')

% Other cumulative values
z(3:2:end,:) = cumsum(a0(1:2:end,:).*y(1:2:m-2,:) +...
    a1(1:2:end,:).*y(2:2:m-1,:) +...
    a2(1:2:end,:).*y(3:2:m,:),1);
z(4:2:end,:) = cumsum(a0(2:2:end,:).*y(2:2:m-2,:) +...
    a1(2:2:end,:).*y(3:2:m-1,:) +...
    a2(2:2:end,:).*y(4:2:m,:),1) +...
    repmat(z(2,:),ceil((m-3)/2),1);

% Resizing
siz = size(y); siz(1) = max(1,siz(1));
z = reshape(z,[ones(1,nshifts),siz]);
if ~isempty(perm), z = ipermute(z,perm); end

end
