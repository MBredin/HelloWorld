function [A,IF,PSI,SIGMA] = IMFdemodVMD(S,varargin)
%==========================================================================
% Call Syntax:  [A,IF,PSI,SIGMA] = IMFdemodVMD(S,fs,L,alpha, tau)
%
% Description:  This function perfoms IMF demodulation.
%
% Input Arguments:
%
%	Name: S
%	Type: matrix (real)
%	Description: matrix of IMF modes  (each as a column), with residual in last column.
%
%   Name: fs
%   Type: scalar
%   Description: sampling freq
%
%	Name: L (optional)
%	Type: integer  [default value:0]
%	Description: ignore L additional point in the vacinity of sigma's zeros
%	crossings to determine the inst. freq, (this help with numerical computational
%	problems)
%
%   Name: alpha
%   Type: 
%   Description: 
%
%	Name: tau
%	Type: 
%	Description: 
%
% Output Arguments:
%
%	Name: A
%	Type: matrix (real)
%	Description: matrix of inst. amplitudes, columns correspond to each IMF
%
%	Name: IF
%	Type: matrix (real)
%	Description: normalized inst freq, columns correspond to each IMF
%
%	Name: PSI
%	Type: matrix (complex)
%	Description: matrix of complex IMF estimate, columns correspond to each IMF
%
%	Name: SIGMA
%	Type: matrix (real)
%	Description: estimated imaginary part of each IMF, columns correspond to each IMF
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
%   [1] "On the HHT, its problems, and some solutions", Reference: Rato, R. T., Ortigueira, M. D., and Batista, A. G.,
%       Mechanical Systems and Signal Processing , vol. 22, no. 6, pp. 1374-1394, August 2008.
%
% Notes:
%
%
% Function Dependencies:    IMFdemod.m    
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

params.IAestMaxIter = 5;
params.IFestThrowOut = 0;
params.fs = 1;
params.VMDdemodAlpha = 5e4;
params.VMDdemodTau = 1e-4;
params.IMFsigmaEstDeriv = 'center3';
params.demodPhaseDeriv = 'center9';

params = parse_pv_pairs(params,varargin);

%-----------
% Initialize
%-----------

K = 1;
DC = false;
init = 1;
tol = 1e-6;

%-----
% Main
%-----

[A,~,PSI,SIGMA] = IMFdemod(S,'fs', params.fs,'IFestThrowOut', params.IFestThrowOut,'IAestMaxIter',params.IAestMaxIter,'IMFsigmaEstDeriv',params.IMFsigmaEstDeriv);
for k = 1:size(S,2);
    [u, ~, ~] = VMD(S(:,k), params.VMDdemodAlpha, params.VMDdemodTau, K, DC, init, tol);
    [~,IF(:,k),~,~] = IMFdemod(transpose(u),'fs', params.fs,'IFestThrowOut', params.IFestThrowOut,'IAestMaxIter',params.IAestMaxIter,'IMFsigmaEstDeriv',params.IMFsigmaEstDeriv);    
end





