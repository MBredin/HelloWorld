function rPMax= parabolicMaxs(aS,opt)
% Call Syntax:  rPMax= rGetPMaxs_s(aS)
%
% Description: This function gets Parabolic Maxs, plateaus out
%
% Input Arguments:
%	Name: aS
%	Type: vector (real)
%	Description: input signal
%
% Output Arguments:
%
%	Name: rPMax
%	Type: vector (real)
%	Description: Parabolic Maxs
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
%
%   This program was derivved from:
%
%       �On the HHT, its problems, and some solutions�, Reference: Rato, R. T., Ortigueira, M. D., and Batista, A. G.,
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

%------------------
% Check valid input
%------------------

if nargin<2
    opt = 'discard';
end






kS= aS(:);
quntLenS=length(kS);
quntMaxCnt=0;


%VECTORIZED BY STEVEN
if (quntLenS>2)     %if signal has enough length
    CNT = 2:(quntLenS-1);
    TrueFalse = ( ((kS(CNT) > kS(CNT+1))) & ((kS(CNT) >= kS(CNT-1))) | ((kS(CNT) >= kS(CNT+1))) & ((kS(CNT) > kS(CNT-1))) );
    kSMNdx1 = find(TrueFalse==true)+1;
    quntMaxCnt = sum(TrueFalse);
    %kSMVal=kS(kSMNdx1);
end


%-----
% Now we have the Maxs, lets get the Parabolic Maxs
intGapMax= max(kS)-min(kS);

JJ = 1:quntMaxCnt;     %for all Maxs

YA = kS(kSMNdx1(JJ)-1);  % Sample point before         (%xa= -1; xb= 0; xc= 1;)
YB = kS(kSMNdx1(JJ));    % Sample point, == kSMVal(jj)
YC = kS(kSMNdx1(JJ)+1);  % Sample point after

D1  = (-4.*YB+2.*YA+2.*YC);
D2 = (-16.*YB+ 8.*YA+ 8.*YC);

XV(D1==0) = kSMNdx1(JJ(D1==0));
XV(D1~=0) = kSMNdx1(JJ(D1~=0))+(YA(JJ(D1~=0))-YC(JJ(D1~=0)))./D1(JJ(D1~=0));

YV(D2(JJ)==0)  = YB(JJ(D2(JJ)==0));
YV(D2(JJ)~=0)  = YB(JJ(D2(JJ)~=0))+ (2*YC(JJ(D2(JJ)~=0)).*YA(JJ(D2(JJ)~=0))- YA(JJ(D2(JJ)~=0)).*YA(JJ(D2(JJ)~=0))- YC(JJ(D2(JJ)~=0)).*YC(JJ(D2(JJ)~=0)))./D2(JJ(D2(JJ)~=0));

if ~isempty(XV)
    TEMP(1) = false;
    TEMP(2:quntMaxCnt) = ( (XV(2:quntMaxCnt)==XV((2:quntMaxCnt)-1))|(abs(YV(2:quntMaxCnt)-YV((2:quntMaxCnt)-1))./abs(XV(2:quntMaxCnt)-XV((2:quntMaxCnt)-1)))> (2*intGapMax) );
    kSPMTim1 = XV;
    kSPMVal  = YV;
    if sum(TEMP)>0
        kSPMTim1(find(TEMP)-1) = (XV(find(TEMP))+ XV(find(TEMP)-1))./2; %#ok<FNDSB>
        kSPMVal(find(TEMP)-1) = max([YV(find(TEMP));YV(find(TEMP)-1)],[],1)'; %#ok<FNDSB>
    end
    kSPMTim1 = kSPMTim1(~TEMP)';
    kSPMVal = kSPMVal(~TEMP)';
else
    kSPMTim1 = [];
    kSPMVal = [];
end

%-----

if quntMaxCnt>0
    if ( kS(1) >= kSPMVal(1) )
        kSPMTim1= [1; kSPMTim1];  kSPMVal=[kS(1); kSPMVal ];    %Start must be included as a Max
    end
    if ( kS(end) >= kSPMVal(end))
        kSPMTim1= [kSPMTim1; quntLenS];  kSPMVal=[kSPMVal; kS(end)];   %End must be included as a Max
    end
end

if quntMaxCnt==0
    if ( kS(1) > kS(2) )
        kSPMTim1= [1; kSPMTim1];  kSPMVal=[kS(1); kSPMVal ];    %Start must be included as a Max
    end
    if ( kS(end) > kS(end-1))
        kSPMTim1= [kSPMTim1; quntLenS];  kSPMVal=[kSPMVal; kS(end)];   %End must be included as a Max
    end
end
if quntMaxCnt<0
    error('rGetPMaxs_s: Invalid MaxCnt value');
end


rPMax= sortrows([kSPMTim1, kSPMVal]);

%REMOVE MINIMA AT START OR END
if strcmp('discard',opt)
    
    if rPMax(1,1)==1
        %warning('ignoring min at start')
        rPMax = rPMax(2:end,:);
    end
    if rPMax(end,1)==length(aS)
        %warning('ignoring min at end')
        rPMax = rPMax(1:end-1,:);
    end
    
end



end

