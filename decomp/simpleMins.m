function sMin = simpleMins(x,opt)
% Call Syntax: sMin = simpleMins(x,opt)
%
% Description: This function finds minima and correspond time indexes
%
% Input Arguments:
%	Name: x
%	Type: vector (real)
%	Description: input signal
%
%	Name: opt
%	Type: string
%	Description: keep or disard minima at signal start and end
%                   'discard' [default]
%                   'keep'
%
% Output Arguments:
%
%	Name: sMax
%	Type: matrix (real)
%	Description: time index in first column, minima values in second column
%
% References:
%
%
% If you use these files please cite the following:
%
%       @article{HSA2015,
%           title={Theory of the Hilbert Spectrum},
%           author={Sandoval, S. and De~Leon, P.~L.~},
%           journal={{Applied and Computational Harmonic Analysis}},
%           year = {\noop{2015}in review},  }
%
%--------------------------------------------------------------------------
% Notes:
%
%--------------------------------------------------------------------------
% Revision History:
%
%==========================================================================

%------------------
% Check valid input
%------------------

if nargin<2
    opt = 'discard';
end


%-----------
% Initialize
%-----------

X = [[NaN;x(1:end-1)],x,[x(2:end);NaN]];
t = (1:length(x))';

%-----
% Main
%-----

%FIND MINIMA TIME INDEX
[~,imin] = min(X,[],2);
imin = (imin==2);

%REMOVE MINIMA AT START OR END
if strcmp('discard',opt)
    
    if imin(1)==true
        %warning('ignoring min at start')
        imin(1) = false;
    end
    if imin(end)==true
        %warning('ignoring min at end')
        imin(end) = false;
    end
    
end

%RESUTLS
sMin = [t(imin),x(imin)];




