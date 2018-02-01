function y = derivApprox(x,fs,method) 
%==========================================================================
% Call Syntax: y = derivApprox(x,fs,method)
%
% Description: This function estimates the approximate computation of a 
%              derivative using various mumerical techniques.  
%
% Input Arguments:
%   Name: x
%   Type: column vector, ot matix of column vectors
%   Description: signal(s) to be differentiated
%
%   Name: fs (optional)
%   Type: scalar 
%   Description: sampling frequency [default = 1] 
%
%   Name: method (optional)
%   Type: string
%   Description: numerical differentiation method: 
%                   'forward'   - use the forward difference
%                   'backward'  - use the backward difference
%                   'center3'   - use the 3-pt stencil central difference 
%                   'center5'   - use the 5-pt stencil central difference 
%                   'center7'   - use the 7-pt stencil central difference 
%                   'center9'   - use the 9-pt stencil central difference 
%                   'center11'  - use the 11-pt stencil central difference 
%                   'center13'  - use the 13-pt stencil central difference 
%                   'center15'  - use the 15-pt stencil central difference [default]
%
% Output Arguments:
%   Name: y
%   Type: column vector or matrix of column vectors
%   Description: Derivativr approximation
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
% References: [1] http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/central-differences/
%             [2] http://web.media.mit.edu/~crtaylor/calculator.html
%
% Notes: See DEMO_derivApprox.m
%
%
% Function Dependencies:  
%
%--------------------------------------------------------------------------
% Author: Steven Sandoval
%--------------------------------------------------------------------------
% Creation Date: May 2017
%
% Revision History: 9-12-2017 S. Terrazas - update function call
%
%==========================================================================

%------------------
% Check valid input
%------------------

if nargin<3
    method = 'center15';
end

if nargin<2
    fs = 1;
end

if not(strcmp(method,'forward')||strcmp(method,'backward')||strcmp(method,'center3')||strcmp(method,'center5')||strcmp(method,'center7')||strcmp(method,'center9')||strcmp(method,'center11')||strcmp(method,'center13')||strcmp(method,'center15')) %check c11 c13 c15
    error('Error (derivApprox): invalid derivative method.')
end

%DEAL WITH MATRIX INPUT
rows=size(x,1);
cols=size(x,2);
if cols>1
    y=zeros(rows, cols);
    for k=1:size(x,2)
        [y(:,k)]=derivApprox(x(:,k), fs,method);
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
    case 'forward'
        y = diff(x);
        y = [NaN;y];
    case 'backward'
        y = diff(x);
        y = [y;NaN];
    case 'center3'
        y =  ( x(3:end)-x(1:end-2) )./2;
        y = [NaN;y;NaN];
    case 'center5'
        y = ( x(1:end-4)-8.*x(2:end-3)+8.*x(4:end-1)-x(5:end) )./12;
        y = [NaN(2,1);y;NaN(2,1)];
    case 'center7'
        y = ( -x(1:end-6)+ 9.*x(2:end-5)- 45.*x(3:end-4) +45.*x(5:end-2) -9.*x(6:end-1) + x(7:end)   )./60;
        y = [NaN(3,1);y;NaN(3,1)];
    case 'center9'
        y = (3.*x(1:end-8) - 32.*x(2:end-7) + 168.*x(3:end-6) -672.*x(4:end-5) + 672.*x(6:end-3) -168.*x(7:end-2) + 32.*x(8:end-1) -3.*x(9:end)  )./840;
        y = [NaN(4,1);y;NaN(4,1)];
    case 'center11'
        y = (-2.*x(1:end-10) +25.*x(2:end-9) -150.*x(3:end-8) +600.*x(4:end-7) -2100.*x(5:end-6)         + 2100.*x(7:end-4)  -600.*x(8:end-3) +150.*x(9:end-2) -25.*x(10:end-1) +2.*x(11:end)  )./2520;
        y = [NaN(5,1);y;NaN(5,1)];
    case 'center13'
        y = (5.*x(1:end-12) - 72.*x(2:end-11) + 495.*x(3:end-10) -2200.*x(4:end-9)   +7425.*x(5:end-8)  -23760.*x(6:end-7)                 +23760.*x(8:end-5)   + -7425.*x(9:end-4)   + 2200.*x(10:end-3) -495.*x(11:end-2) + 72.*x(12:end-1) -5.*x(13:end)  )./27720;
        y = [NaN(6,1);y;NaN(6,1)];
    case 'center15'
        y = (-67553994410557440.*x(1:end-14) +1103381908705771500.*x(2:end-13) -8606378887905019000.*x(3:end-12) +43031894439525120000.*x(4:end-11) -157783612944925330000.*x(5:end-10)  + 473350838834776000000.*x(6:end-9)  -1.4200525165043282e21.*x(7:end-8)             +1.4200525165043282e21.*x(9:end-6)-473350838834776000000.*x(10:end-5)+157783612944925330000.*x(11:end-4)-43031894439525120000.*x(12:end-3)  +8606378887905019000.*x(13:end-2) -1103381908705771500.*x(14:end-1) +67553994410557440.*x(15:end)  )./1.622917161719232e21;
        y = [NaN(7,1);y;NaN(7,1)];
end
y = y.*fs;

end
