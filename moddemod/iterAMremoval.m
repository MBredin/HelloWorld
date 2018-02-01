function [A,S_FM,A1] = iterAMremoval(S,K)
%==========================================================================
% Call Syntax:  [A,S_FM,A1] = iterAMremoval(S,K)
%
% Description:  
%
% Input Arguments:
%	Name: S
%	Type: matrix (real)
%	Description: matrix of IMF modes  (each as a column), with residual in last column.
%
%	Name: K
%	Type: integer (optional)
%	Description: number of iterations [default value: 3]
%
% Output Arguments:
%
%	Name: A
%	Type: matrix (real)
%	Description: matrix of inst. amplitudes, after K iterations.
%
%	Name: S_FM
%	Type: matrix (real)
%	Description: matrix of amplitude normalized IMF modes, columns correspond to each IMF
%
%	Name: A1
%	Type: matrix (real)
%	Description: matrix of inst. amplitudes, after 1 iteration. You may
%	want to use this value as the IA esitimate.
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
%   [1] R. Rato, M. Ortigueira, and A. Batista, “On the HHT, its problems, and
%       some solutions," Mechanical Systems and Signal Processing, vol. 22,
%       no. 6, pp. 1374–1394, 2008.
%
%   [2] Huang, Norden E., et al. "On instantaneous frequency." Advances in
%       Adaptive Data Analysis 1.02 (2009): 177-229.
%
% Notes:
%
%
% Function Dependencies:    IAest.m
%
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

if nargin<2; 
    K=5; 
end

%-----------
% Initialize
%-----------

A = ones(size(S));
S_FM = S;

%-----
% Main
%-----

for k = 1:K                 %loop over iterations
    Atmp = IAest(S_FM);     %perform IA estimation
    if k==1                 %if first iteration
        A1 = Atmp;          %save IA
    end                     %end if first iteration
    A = A.*Atmp;            %compute IA estimation oat the k-th iteration
    S_FM = S_FM./Atmp;      %compute the amplitude normalized signal at the k-th iteration
end                         %endloop over iterations


