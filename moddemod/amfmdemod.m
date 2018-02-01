function [A,IF,S,SIGMA] = amfmdemod(PSI,fs,diffMode)
%==========================================================================
% Call Syntax: [A,IF,S,SIGMA] = amfmdemod(PSI,fs,diffMode)
%
% Description:  This function demodulates a complex AM-FM component.
%
% Input Arguments:
%   Name: PSI
%   Type: complex matrix (or vector)
%   Description: each column is a complex AM-FM component
%
%   Name: fs
%   Type: scalar
%   Description: sampling freq
%
%   Name: diffMode
%   Type: string
%   Description: numerical differentiation method: 'forward' or 'backward' or
%   'center3' or 'center5' or 'center7' or 'center9'[default]
%
% Output Arguments:
%   Name: A
%   Type: real matrix (or vector)
%   Description: each column is the instantaneous amplitude
%
%   Name: IF
%   Type: real matrix (or vector)
%   Description: each column is the instantaneous frequency
%
%   Name: S
%   Type: real matrix (or vector)
%   Description: each column is the real part of the AM-FM signal
%
%   Name: SIGMA
%   Type: real matrix (or vector)
%   Description: each column is the imaginary part of the AM-FM signal
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
% Function Dependencies:    derivApprox.m
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

if nargin<3
    diffMode  = 'center9';
end


%-----------
% Initialize
%-----------

IF = zeros(size(PSI)); %allocate memeory space


%-----
% Main
%-----

for j= 1:size(PSI,2);                               %loop over each column in PSI
    thetaU = unwrap(angle(PSI(:,j)));               %unwrap the phase function
    fi = derivApprox(thetaU,fs,diffMode)./(2*pi);   %compute the instantaneous frequency in Hz
    IF(:,j) = fi(:);                                %store result as a column vector 
end                                                 %end loop over each column in PSI

A = abs(PSI);       %extract the instantaneous amplitude
S = real(PSI);      %extract the real signal part
SIGMA = imag(PSI);  %extract the imaginary signal part





