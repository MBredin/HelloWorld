function rPMinExtrapol = rPMinExtrapol_s(rPM, rPm, quntL)
% Call Syntax:  rPMinExtrapol = rPMinExtrapol_s(rPM, rPm, quntL)
%
% Description: This function performs Time-mirrored bottom extrema (Parabolic Mins) extrapolation
%
% Input Arguments:
%	Name: rPM
%	Type:
%	Description:
%
%	Name: rPm
%	Type:
%	Description:
%
%	Name: quntL
%	Type:
%	Description:
%
% Output Arguments:
%
%	Name: rPMaxExtrapol
%	Type: vector (real)
%	Description: Time-mirrored top extrema (Parabolic Mins) extrapolation
%
% References:
%
% If you use these files please cite the following:
%
%       @article{HSA2015,
%           title={Theory of the Hilbert Spectrum},
%           author={Sandoval, S. and De~Leon, P.~L.~},
%           journal={{Applied and Computational Harmonic Analysis}},
%           year = {\noop{2015}in review},  }
%
%
%   This program was derivved from:
%
%       “On the HHT, its problems, and some solutions”, Reference: Rato, R. T., Ortigueira, M. D., and Batista, A. G.,
%       Mechanical Systems and Signal Processing , vol. 22, no. 6, pp. 1374-1394, August 2008.
%       Authors: Raul Rato (rtr@uninova.DOT.pt) and Manuel Ortigueira (mdortigueira@uninova.pt or mdo@fct.unl.pt)
%

%--------------------------------------------------------------------------
% Notes:
%
%--------------------------------------------------------------------------
% Revision History:
%
%   History:    V1.00 First version (R. Rato)
%               V1.01 Count mismatch detection increased from 1 to 2 (R. Rato)
%               V2.00 vectorize (S. Sandoval)
%               V3.00 myHHT(S. Sandoval)
%
% WARNING: This software is a result of our research work and is supplied without any garanties.
%           We would like to receive comments on the results and report on bugs.
%
%==========================================================================


    %------------------------------
    %Init
    %------------------------------
    rPM= sortrows(rPM); %assumes nothing on rPM sort order
    rPm= sortrows(rPm); %assumes nothing on rPm sort order
    
    kTopTim1= rPM(:,1); kTopVal= rPM(:,2);
    kDwnTim1= rPm(:,1); kDwnVal= rPm(:,2);
    
    %------------------------------
    %Start extrapolation
    %------------------------------
    if ( (kTopTim1(1)== 1) && (kDwnTim1(1)== 1) )
        %disp ('            rPMinExtrapol: Poliextrema at signal''s start');
    elseif ( (kTopTim1(1)<1) || (kDwnTim1(1)< 1) )
        %disp ('            rPMinExtrapol: Invalid extrema at signal''s start');
    else
        kDwnTim1=[2-kTopTim1(1); kDwnTim1];     % New first Dwn at the (one based) specular Max
        kDwnVal=[kDwnVal(1); kDwnVal];          % Same Val as old first Dwn
    end
    
    %------------------------------
    % End extrapolation
    %------------------------------
    if ( (kTopTim1(end)== quntL) && (kDwnTim1(end)== quntL) )
        %disp ('            rPMinExtrapol: Poliextrema at signal''s end');
    elseif ( (kTopTim1(end)> quntL) || (kDwnTim1(end)> quntL) )
        %disp ('            rPMinExtrapol: Invalid extrema at signal''s end');
    else
        kDwnTim1=[kDwnTim1; (2*quntL - kTopTim1(end))];     % New last Dwn at the specular Max
        kDwnVal=[ kDwnVal; kDwnVal(end)];          % Same Val as old last Dwn
    end
    
    %------------------------------
    % return value
    %------------------------------
    rPMinExtrapol= sortrows([kDwnTim1, kDwnVal]);
    
end

