function [varphi,numIter,WBiasAVGdB] = siftNEW(r,varargin)
% Call Syntax:  
%
% Description: This function perfoms sifting to remove a single IMF
%
% Input Arguments:
%	Name: r
%	Type: vector (real)
%	Description: original signal from which to extract one IMF estimate
%
%	Name: StopThresh
%	Type: scalar (positive)
%	Description: Sifting stop criterion 
%
%	Name: alpha
%	Type: scalar (positive)
%	Description: step size
%
%	Name: maxSifts (optional)
%	Type: integer (positive)
%	Description: Sifting stop criterion [default=50]
%
% Output Arguments:
%
%	Name: varphi
%	Type: vector (real)
%	Description: vector with the first IMF
%
%==========================================================================

%------------------
% Check valid input
%------------------

%DEFAULT PARAMETERS
params.SiftStopThresh = -120;
params.maxSifts = 500;
params.siftAlpha = 1;
params.VerboseSift = false;
params.ExtremaID = 'parabolic';
params.postextrapolate = 'none';
params.tangentInterpIter = 0;
params.fs = 1;
params.extremaEndpoints = 'discard';
params.lpcExtrapLen = 50;
params.numLPC = 200;
params.tangentDeriv = 'center9'; 
params.siftInterp = 'spline';
params.zeroMeanIMFs = true;

%GET USER SPECIFED PARAMETETERS
params = parse_pv_pairs(params,varargin);


%-----------
% Initialize
%-----------
if (params.lpcExtrapLen>0)
    varphi = [lpredict2(r, params.numLPC, params.lpcExtrapLen, 'pre');r;lpredict2(r, params.numLPC, params.lpcExtrapLen, 'post')];
else
    varphi = r;                      %IMF estimate starts at the input signal
end
IMF_len = length(varphi);            %get signal length
WBiasAVGdB = NaN(params.maxSifts,1); %initialize vector for convergence measure
e = zeros(length(varphi),1);         %initialize the trend signal
numIter = 0;                         %initialize sifting iteration counter

%-----
% Main
%-----

while (true)
    numIter = numIter+1;
    
    %SIFTING
    varphi = varphi- params.siftAlpha*e;   

    %subtract alpha*bias from Imf estimate
    [u_pAndt_p, l_qAndt_q] = getExtrama(varphi,'method', 'parabolic','postextrapolate', 'oddReflect','extremaEndpoints', 'discard');    %get maxima and minima 
    
    for i_loop = 1:params.tangentInterpIter
        [u_pAndt_p,l_qAndt_q] = tangentInterPts(varphi,params.fs,u_pAndt_p,l_qAndt_q,params.tangentDeriv);
    end
    if strcmp(params.siftInterp,'spline')
        u = spline(u_pAndt_p(:,1), u_pAndt_p(:,2), 1:IMF_len);                                  %interpolate maxima to get upper envelope
        l = spline(l_qAndt_q(:,1), l_qAndt_q(:,2), 1:IMF_len);                                  %interpolate minima to get lower envelope
    else
        error('sift: invalid sifting iterpolator');
    end
    e = transpose((u+l)/2);                                                                 %average the upper and lower envelopes to get the bias
    
    %CHECK FOR STOP CONDITIONS
    WBias= transpose(e(params.lpcExtrapLen+1:end-params.lpcExtrapLen) )*transpose(e(params.lpcExtrapLen+1:end-params.lpcExtrapLen) )';
    WBiasAVGdB(numIter) = 10*log10(WBias/length(varphi));                                                                                           %compute average bias enegy in dB
    if (params.VerboseSift)                                                                                                                         %if verbose
    disp(['        Iteration: ',num2str(numIter),', Average Bias Energy (dB): ',num2str(WBiasAVGdB(numIter)),', Target: ',num2str(params.SiftStopThresh)]) %dispay progress
    end                                                                                                                                             %end if verbose
    if (WBiasAVGdB(numIter)<params.SiftStopThresh)                                                                                                  %if stop threshold is reached
        break;                                                                                                                                      %then break
    end                                                                                                                                             %end if stop threshold is reached
    if (numIter>=params.maxSifts)                                                                                                                   %if max iterations reached
        varphi = varphi- params.siftAlpha*e;                                                                                                            %subtract alpha*bias from Imf estimate
        if (params.VerboseSift)
            warning('sift: forced exit, max sifting iterations reached.');                                                                            %issue warning
        end
        break;                                                                                                                                      %break
    end                                                                                                                                             %end if max iterations reached
    
end
WBiasAVGdB = WBiasAVGdB(1:numIter); %return bias energy only for iterations that were run
if (params.lpcExtrapLen>0)
    varphi = varphi(params.lpcExtrapLen+1:end-params.lpcExtrapLen);
end
if (params.zeroMeanIMFs)
    varphi = varphi-mean(varphi);
end


