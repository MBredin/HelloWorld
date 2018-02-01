function [A,IF,PSI,SIGMA] = IMFdemod(S,varargin)
%==========================================================================
% Call Syntax:  [A,IF,PSI,SIGMA] = IMFdemod(S,fs,L)
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
% Function Dependencies:    iterAMremoval.m
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
params.IMFsigmaEstDeriv = 'center3';
params.demodPhaseDeriv = 'center9';

params = parse_pv_pairs(params,varargin);

%-----------
% Initialize
%-----------

S_FM    = zeros(size(S)); %allocate memory
A       = zeros(size(S)); %allocate memory
SIGMA_FM   = zeros(size(S)); %allocate memory
SGN     = zeros(size(S)); %allocate memory


%-----
% Main
%-----

for j = 1:size(S,2)
    
    [~,S_FM(:,j),A(:,j)] = iterAMremoval(S(:,j),params.IAestMaxIter);               %compute real part of the amplitude normailed (FM) signal
    SIGMA_FM(:,j)        = real( sqrt(1-S_FM(:,j).^2));         %estimate the magnitude of the imaginary signal part
    SGN(:,j)         = -sign(  derivApprox(S_FM(:,j),params.fs,params.IMFsigmaEstDeriv)     );              %take derivative of the real signal part
    SIGMA_FM(:,j)        = SGN(:,j).*SIGMA_FM(:,j);                %append sign to the estimate of the imaginary signal part
    
    nanLoc = crossing(SIGMA_FM(:,j));                              %find the zero crossing of the imaginary signal part
    SIGMA_FM(nanLoc,j)=NaN;                                        %discard samples at the zero crossing
    if params.IFestThrowOut>0                                                      %if number of samples to discard is greater than zero (adjacent to the crossing)
        for k = 1:params.IFestThrowOut                                             %loop over number of sample to discard from each side of the crossing
            nanLocPlus = (nanLoc+k);                            %compute indices k samples to the left of the crossing
            nanLocPlus = nanLocPlus( nanLocPlus<size(SIGMA_FM,2) );%ensure indicies are not out of range
            SIGMA_FM(nanLocPlus,j)=NaN;                            %discard valid indicies to the left
            nanLocMinus = (nanLoc-k);                           %compute indices k samples to the right of the crossing
            nanLocMinus  = nanLocMinus (nanLocMinus >0);        %ensure indicies are not out of range
            SIGMA_FM(nanLocMinus,j)=NaN;                           %discard valid indicies to the left
        end                                                     %loop over number of sample to discard from each side of the crossing
    end                                                         %end if number of samples to discard is greater than zero (adjacent to the crossing)
    
    YY = spline(   find(  not(isnan(SIGMA_FM(:,j)))  )  , SIGMA_FM(not(isnan(SIGMA_FM(:,j))),j),   find(isnan( SIGMA_FM(:,j) )) ); %replace dicarded values with interpolated values

    SIGMA_FM(isnan(SIGMA_FM(:,j)),j) = YY;        %save estimated imaginary signal part
    PSI_FM(:,j) = S_FM(:,j) +1i.*SIGMA_FM(:,j);   %save estimated complex value
    [~,IF(:,j)] = amfmdemod(PSI_FM(:,j),'fs',params.fs,'demodPhaseDeriv',params.demodPhaseDeriv);    %compute the estimated instantaneous frequency
    
    
    
    %figure()
    %hold on;
    %plot(S_FM)
    %plot(SIGMA,'k')
    %plot(nanLoc,SIGMA(nanLoc),'mx')
    %temp = NaN(size(SIGMA(:,j)));
    %temp(isnan(SIGMA(:,j)))=YY;
    %plot(temp,'m')
    %plot(IF*100,'r')

end
SIGMA = SIGMA_FM.*A;
PSI = PSI_FM.*A;








