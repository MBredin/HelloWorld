function [S,A,IF,PSI,SIGMA] = HSA_EMDnew(x,varargin)

%------------------
% Check valid input
%------------------ 

%GET PRECONFIG
for i=1:length(varargin)/2
    if strcmpi('preConfig',varargin{2*i-1})
        params.preConfig = varargin{2*i};
    else
        params.preConfig = 'default';
    end
end

%GET USER SPECIFED PARAMETETERS
params = parse_pv_pairs(preConfig(params.preConfig),varargin);


%-----------
% Initialize
%-----------

%REMOVE MEAN
x_mu = mean(x);
x = x-x_mu;

%NORMALIZE
if (params.computeNormalized)
    x_max = max(abs(x));
    x = x./x_max;
end

%INITIALIZE VARIABLES
r = x(:);
k = 0;
Wx = r'*r;                       %original signal energy
S = [];
PREV_cutoff = 0.5;
b = [];
a = [];

%OPEN PARPOOL
if (params.openParPool)&&(params.trials>1)
    if isempty(gcp('nocreate'))
        myCluster  = parcluster();                           %Get the current parallel cluster information
        myPool = parpool(myCluster,myCluster.NumWorkers);    %Open a parallel pool with as many workers are available on the cluster
        infoStr = ['NumWorkers: ',num2str(myPool.NumWorkers)];
        disp(infoStr);
    end
end

%RESET RANDOM 
if (params.resetSeed)
   rng('default');
end

warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary');
%-----
% Main
%-----

if (params.statusEMD);disp('Decomposing Signal...');end %disp progress
while (true)
    
    
    k = k+1;
    
    if (params.statusEMD);disp(['    Computing IMF ',num2str(k),'...']); end    %disp progress
    
    
    %ESTIMATE THE kth IMF FROM RESIDUE
    
    if (params.trials == 1)
        %EMD
        [varphi,~,~] = siftNEW(r,   'SiftStopThresh', params.SiftStopThresh, 'maxSifts', params.maxSifts, 'siftAlpha', params.siftAlpha, ...
            'VerboseSift', params.VerboseSift, 'ExtremaID', params.ExtremaID, 'postextrapolate', params.postextrapolate, ...
            'tangentInterpIter', params.tangentInterpIter, 'fs', params.fs, 'extremaEndpoints', params.extremaEndpoints, ...
            'lpcExtrapLen', params.lpcExtrapLen, 'numLPC', params.numLPC, 'tangentDeriv', params.tangentDeriv, 'siftInterp', params.siftInterp, ...
            'zeroMeanIMFs', params.zeroMeanIMFs);
    
    else
        
        
        if strcmp(params.maskSignal,'filteredNoise')&& (k>1)
            %DEMODULATE THE IMFs
            if strcmp(params.demodMethod,'IMFdemod')
                [a_k,if_k] = IMFdemod(S(:,k-1),'fs', params.fs,'IFestThrowOut', params.IFestThrowOut,'IAestMaxIter',params.IAestMaxIter,...
                    'IMFsigmaEstDeriv',params.IMFsigmaEstDeriv,'demodPhaseDeriv',params.demodPhaseDeriv);
            elseif strcmp(params.demodMethod,'VMDdemod')
                [a_k,if_k] = IMFdemodVMD(S(:,k-1),'fs', params.fs,'IFestThrowOut', params.IFestThrowOut,'IAestMaxIter',params.IAestMaxIter, ...
                    'VMDdemodAlpha',params.VMDdemodAlpha,'VMDdemodTau',params.VMDdemodTau,...
                    'IMFsigmaEstDeriv',params.IMFsigmaEstDeriv,'demodPhaseDeriv',params.demodPhaseDeriv);
            elseif strcmp(params.demodMethod,'HTdemod')
                [a_k,if_k] = HTdemod(S(:,k-1),'demodPhaseDeriv',params.demodPhaseDeriv,'fs', params.fs);
            end
            
            %FILTER DESIGN
            if strcmp(params.maskSignalcutoff,'IAweightedIF')                                                    %AMPLITUDE WEIGHTED MEAN INST FREQ
                cutoff = nanmean(sum(if_k.*a_k) ./ sum(a_k))/params.fs;                                                            %amplitude weighted mean IF
                if not(isfinite(cutoff) && (cutoff>0) && (cutoff<0.5))                                %check for valid cutoff freq
                    cutoff = 0.75* PREV_cutoff;
                    if (params.VerboseEMD); warning('       HSA: Mean IF undefined for previous IMF'); end            %if bad, give warning and correct
                end                                                                                                 %end if bad
                [b,a] = butter(params.maskFilterOrder,params.maskFreqSpacing*cutoff);                                                     %filter noise
                PREV_cutoff = cutoff;                                                                             %save cutoff freq
            elseif strcmp(params.maskSignalcutoff,'EweightedIF')                                                   %ENERGY WEIGHTED MEAN INST FREQ
                cutoff = nanmean(sum(if_k.^2.*a_k) ./ sum(if_k.*a_k))/params.fs;                                               %energy weighted mean IF
                if not(isfinite(cutoff) && (cutoff>0) && (cutoff<0.5))                                           %check for valid cutoff freq
                    cutoff = 0.75* PREV_cutoff;
                    if (params.VerboseEMD);warning('       HSA: Mean IF undefined for previous IMF');end              %if bad, give warning and correct
                end                                                                                                 %end if bad
                [b,a] = butter(params.maskFilterOrder,params.maskFreqSpacing*cutoff);                                                     %filter noise
                PREV_cutoff = cutoff;                                                                             %save cutoff freq
            elseif strcmp(params.maskSignalcutoff,'meanIF')                                                       %MEAN INST FREQ
                cutoff = nanmean(if_k)/params.fs;                                                                            %mean IF
                if not(isfinite(cutoff) && (cutoff>0) && (cutoff<0.5))                                           %check for valid cutoff freq
                    cutoff = 0.75* PREV_cutoff;
                    if (params.VerboseEMD);warning('       HSA: Mean IF undefined for previous IMF'); end            %if bad, give warning and correct
                end                                                                                                 %end if bad
                [b,a] = butter(params.maskFilterOrder,params.maskFreqSpacing*cutoff);                                                     %filter noise
                PREV_cutoff = cutoff;                                                                             %save cutoff freq
            end
        end
        
        
        
        %CEEMD w/TM
        parfor trialCounter = 1:params.trials
            
            %DESIGN MASKING SIGNAL
            warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary')
            
            
            if strcmp(params.maskSignal,'tone') %
                maskSignal = params.maskAmplitude .* sin(pi*params.maskFreqSpacing^k*(1:length(r))'+ 2*pi*rand(1) );
                
                
            elseif strcmp(params.maskSignal,'siftedNoise')
                if strcmp(params.maskNoiseType,'uniform')
                    maskSignal = params.maskAmplitude * (2.*rand(size(r))-1 ) ;                          %filter noise
                elseif strcmp(params.maskNoiseType,'normal')
                    maskSignal = params.maskAmplitude * randn(size(r));                          %filter noise
                end
                
                for m = 1:k-1
                    maskSignal = siftNEW(maskSignal,   'SiftStopThresh', params.SiftStopThresh, 'maxSifts', params.maxSifts, 'siftAlpha', params.siftAlpha, ...
                        'VerboseSift', params.VerboseSift, 'ExtremaID', params.ExtremaID, 'postextrapolate', params.postextrapolate, ...
                        'tangentInterpIter', params.tangentInterpIter, 'fs', params.fs, 'extremaEndpoints', params.extremaEndpoints, ...
                        'lpcExtrapLen', params.lpcExtrapLen, 'numLPC', params.numLPC, 'tangentDeriv', params.tangentDeriv, 'siftInterp', params.siftInterp, ...
                        'zeroMeanIMFs', params.zeroMeanIMFs)
                end
                
                
            elseif strcmp(params.maskSignal,'filteredNoise')
                
                if k == 1                                                                           %if first iteration
                    if strcmp(params.maskNoiseType,'uniform')
                        maskSignal =  params.maskAmplitude * (2.*rand(size(r))-1 )  ;                          %filter noise
                    elseif strcmp(params.maskNoiseType,'normal')
                        maskSignal = params.maskAmplitude * randn(size(r))   ;                          %filter noise
                    end
                    
                else                                                                                %else
                    
                    
                    if strcmp(params.maskNoiseType,'uniform')
                        maskSignal = filter(b,a,  params.maskAmplitude * (2.*rand(size(r))-1 )  );                          %filter noise
                    elseif strcmp(params.maskNoiseType,'normal')
                        maskSignal = filter(b,a,  params.maskAmplitude * randn(size(r))   );                          %filter noise
                    end
                    
                end                                                                                 %endif
                
            end
            
            
            
            [varphiPlus,~,~]  = siftNEW(r+maskSignal,   'SiftStopThresh', params.SiftStopThresh, 'maxSifts', params.maxSifts, 'siftAlpha', params.siftAlpha, ...
                'VerboseSift', params.VerboseSift, 'ExtremaID', params.ExtremaID, 'postextrapolate', params.postextrapolate, ...
                'tangentInterpIter', params.tangentInterpIter, 'fs', params.fs, 'extremaEndpoints', params.extremaEndpoints, ...
                'lpcExtrapLen', params.lpcExtrapLen, 'numLPC', params.numLPC, 'tangentDeriv', params.tangentDeriv, 'siftInterp', params.siftInterp, ...
                'zeroMeanIMFs', params.zeroMeanIMFs);
        
            [varphiMinus,~,~] = siftNEW(r-maskSignal,   'SiftStopThresh', params.SiftStopThresh, 'maxSifts', params.maxSifts, 'siftAlpha', params.siftAlpha, ...
                'VerboseSift', params.VerboseSift, 'ExtremaID', params.ExtremaID, 'postextrapolate', params.postextrapolate, ...
                'tangentInterpIter', params.tangentInterpIter, 'fs', params.fs, 'extremaEndpoints', params.extremaEndpoints, ...
                'lpcExtrapLen', params.lpcExtrapLen, 'numLPC', params.numLPC, 'tangentDeriv', params.tangentDeriv, 'siftInterp', params.siftInterp, ...
                'zeroMeanIMFs', params.zeroMeanIMFs);
            
            temp = 1/2 .* (varphiPlus + varphiMinus);
            
            %CHECK FOR PROBLEMS
            if (sum(abs(temp)>1e3)>0)||(sum(not(isfinite(temp)))>0)     %look for huge values or non-finite values in solution
                tempIMFs(:,trialCounter) = NaN(size(temp));                       %if present, throw out trail
            else                                                        %else
                tempIMFs(:,trialCounter) = temp;                                  %proceed
            end
            
        end
        temp = nanmean(tempIMFs,2); %average the ensemble
        varphi = siftNEW(temp,   'SiftStopThresh', params.SiftStopThresh, 'maxSifts', params.maxSifts, 'siftAlpha', params.siftAlpha, ...
                    'VerboseSift', params.VerboseSift, 'ExtremaID', params.ExtremaID, 'postextrapolate', params.postextrapolate, ...
                    'tangentInterpIter', params.tangentInterpIter, 'fs', params.fs, 'extremaEndpoints', params.extremaEndpoints, ...
                    'lpcExtrapLen', params.lpcExtrapLen, 'numLPC', params.numLPC, 'tangentDeriv', params.tangentDeriv, 'siftInterp', params.siftInterp, ...
                    'zeroMeanIMFs', params.zeroMeanIMFs);
        
    end
    S = [S,varphi(:)];
    if (params.statusEMD); disp(['      IMF ',num2str(k),' complete.']); end %disp progress
    
    %UPDATE RESIDUE
    r = r - varphi(:);
    
    %CHECK STOP CONDITIONS
    numOscil = countOscil(r);               %get oscilation count
    if (r'*r)>0                             %if
        ThreshCheck = 10*log10(Wx/(r'*r));  %residual energy
    else                                    %otherwise
        ThreshCheck = Inf;                  %residual energy
    end                                                                                                                                            %end if stop threshold is reached
    if ((ThreshCheck>=params.EMDStopThresh) || (numOscil<=2) )
        if ((r'*r))>(10^-12)
            S = [S,r];
        end
        break;
    end
    
end

%UNNORMALIZE
if (params.computeNormalized)
    S = S.*x_max;
end

%DEMODULATE THE IMFs
if strcmp(params.demodMethod,'IMFdemod')
    [A,IF,PSI,SIGMA] = IMFdemod(S,'fs', params.fs,'IFestThrowOut', params.IFestThrowOut,'IAestMaxIter',params.IAestMaxIter,...
                                'IMFsigmaEstDeriv',params.IMFsigmaEstDeriv,'demodPhaseDeriv',params.demodPhaseDeriv);                         
elseif strcmp(params.demodMethod,'VMDdemod') 
    [A,IF,PSI,SIGMA] = IMFdemodVMD(S,'fs', params.fs,'IFestThrowOut', params.IFestThrowOut,'IAestMaxIter',params.IAestMaxIter, ...
                                    'VMDdemodAlpha',params.VMDdemodAlpha,'VMDdemodTau',params.VMDdemodTau,...
                                    'IMFsigmaEstDeriv',params.IMFsigmaEstDeriv,'demodPhaseDeriv',params.demodPhaseDeriv);                              
elseif strcmp(params.demodMethod,'HTdemod')   
    [A,IF,PSI,SIGMA] = HTdemod(S,'demodPhaseDeriv',params.demodPhaseDeriv,'fs', params.fs);
end

%FILTER THE IFs
if (params.demodFiltFreq~=Inf)
    [b,a] = butter(6,params.demodFiltFreq/(1/2));
    for k = 1:size(IF,2)
        IFtemp = IF(:,k);
        IFtemp(isnan(IFtemp)) = nanmean(IFtemp);
        IFtemp = filtfilt(b, a, IFtemp);
        IF(:,k) = IFtemp(:);
    end
end

%RESTORE THE DC COMPONENT
if (params.keepDC)
    S       =  [S,x_mu.*ones(size(x(:)))];
    A       =  [A,x_mu.*ones(size(x(:)))];
    IF      =  [IF,zeros(size(x(:)))];
    PSI     =  [PSI,x_mu.*ones(size(x(:)))];
    SIGMA   =  [SIGMA,zeros(size(x(:)))];
end

%THROW OUT LOW ENERGY IMFs
if (params.EthreshPCT<100)
    ncompStart = size(A,2);
    [Idx, pctOut]  = IAthresh(A,params.EthreshPCT);
    S     = S(:,Idx);
    A     = A(:,Idx);
    IF    = IF(:,Idx);
    PSI   = PSI(:,Idx);
    SIGMA = SIGMA(:,Idx);
    disp(['    ',num2str(size(A,2)),' of ',num2str(ncompStart),' kept, constituting ',num2str(pctOut,15),'% of energy in component set.'])
    
end


