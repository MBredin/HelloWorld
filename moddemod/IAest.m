function A = IAest(S,method)
%==========================================================================
% Call Syntax:  A = IAest(S,method)
%
% Description:  This function perfoms direct amplitude esitmation using the same
%               assumptions made when performing EMD.
%
% Input Arguments:
%	Name: S
%	Type: matrix (real)
%	Description: matrix of IMF modes  (each as a column), with residual in last column.
%
%   Name: method
%	Type: string
%	Description: interpolation method
%           'pchip'  - use Piecewise Cubic Hermite Interpolating Polynomial Interpolation [default]
%           'spline' - use Cubic Spline Interpolation
%
% Output Arguments:
%
%	Name: A
%	Type: matrix (real)
%	Description: matrix of instantaneous amplitude estimations (each as a column)
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
% Notes:
%
%
% Function Dependencies:    parabolicMaxs.m 
%                           parabolicMins.m 
%                           extrapolMaxs.m
%
%--------------------------------------------------------------------------
% Author: Steven Sandoval
%--------------------------------------------------------------------------
% Creation Date: July 2017
%
% Revision History:
%
%==========================================================================

%------------------
% Check valid input
%------------------

if nargin<2
    method = 'pchip';
end

if not(strcmp(method,'pchip')||strcmp(method,'spline'))
    error('Error (IAest): invalid interpolation method.')
end


%-----------
% Initialize
%-----------

A = zeros(size(S)); %allocate memory

%-----
% Main
%-----

for k = 1:size(S,2) %loop over each column  of S

    g = abs(S(:,k));                                %take absolute value of the real signal part
    rPMOri = parabolicMaxs(g);                      %compute maxima of the absolute value of the real signal part
    rPmOri = parabolicMins(g);                      %compute minima of the absolute value of the real signal part
    M = extrapolMaxs(rPMOri, rPmOri, length(g));    %extrapolate maxima of the absolute value of the real signal part
    
    if strcmp(method,'pchip')
        A(:,k) = pchip(M(:,1),M(:,2),1:length(g));  %interpolate using Piecewise Cubic Hermite Interpolating Polynomial
    elseif strcmp(method,'spline')
        A(:,k) = spline(M(:,1),M(:,2),1:length(g)); %interpolate using Cubic Spline
    end
    
end


