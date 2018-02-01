function sMax = simpleMaxs(x,opt)
% Call Syntax: sMax = simpleMaxs(x)
%
% Description: This function finds mixima and correspond time indexes
%
% Input Arguments:
%	Name: x
%	Type: vector (real)
%	Description: input signal
%
%	Name: opt
%	Type: string
%	Description: keep or disard maxima at signal start and end   
%                   'discard' [default]
%                   'keep'
%
% Output Arguments:
%
%	Name: sMax
%	Type: matrix (real)
%	Description: time index in first column, maxima values in second column
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

%FIND MAXIMA TIME INDEX
[~,imax] = max(X,[],2);
imax = (imax==2);

%REMOVE MAXIMA AT START OR END
if strcmp('discard',opt)
    
    if imax(1)==true
        %warning('ignoring max at start')
        imax(1) = false;
    end
    if imax(end)==true
        %warning('ignoring max at end')
        imax(end) = false;
    end
    
end

%RESUTLS
sMax = [t(imax),x(imax)];




