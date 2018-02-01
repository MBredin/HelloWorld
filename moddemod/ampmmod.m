function [psi,s,sigma,fi,t] = ampmmod(a,M,fc,fs,phi,method)
%==========================================================================
% Call Syntax: [psi,s,sigma,fi,t] = ampmmod(a,M,fc,fs,phi,method)
%
% Description:  This function synthesizes an AM-FM component based on a
%               specified instantaneous amplitude, PM message, frequency 
%               reference, sampling frequency, and phase reference.
%
% Input Arguments:
%   Name: a
%   Type: vector
%   Description: AM message
%
%   Name: M
%   Type: vector
%   Description: PM message
%
%   Name: fc
%   Type: scalar
%   Description: center freq (carrier)
%
%   Name: fs
%   Type: scalar
%   Description: sampling freq
%
%   Name: phi
%   Type: scalar
%   Description: initial phase
%
%   Name: method
%   Type: string
%   Description: numerical differentiation method: 'forward' or 'backward' or
%   'center3' or 'center5' or 'center7' or 'center9'[default]
%
% Output Arguments:
%   Name: psi
%   Type: vector (complex)
%   Description: AM-FM component
%
%   Name: s
%   Type: vector (real)
%   Description: real part of the AM-FM signal
%
%   Name: sigma
%   Type: vector (real)
%   Description: imaginary part of the AM-FM signal
%
%   Name: fi
%   Type: vector (real)
%   Description: instantaneous frequency of the AM-FM signal
%
%   Name: t
%   Type: vector (real)
%   Description: time index
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
% Function Dependencies:  derivApprox.m
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

a = a(:); %force column vector
M = M(:); %force column vector

if (nargin ~= 4)&&(nargin ~= 5)&&(nargin ~= 6)
    error('Error (amfmmod): must have 4 to 6 input arguments.');
end;

if (nargin <5)
    phi = 0;
end;

if nargin<6
    method = 'center9';
end

if length(a)~=length(M)
    error('Error (ampmmod): a and M must be of same length')
end


%-----------
% Initialize
%-----------
t = (0:length(m)-1)./fs;    %make time index
t = t(:);                   %force column vector


%-----
% Main
%-----

%SYNTHESIZE AM-FM COMPONENT
theta =  2*pi.*fc.*t + 2*pi.*M + phi;   %compute phase function
psi = a.*exp(1i*theta);                 %compute complex AM-FM component
s = real(psi);                          %extract real part of complex AM-FM component
sigma = imag(psi);                      %extract imaginary part of complex AM-FM component
fi = fc + derivApprox(M,fs,method);     %compute the instananeous frequency


