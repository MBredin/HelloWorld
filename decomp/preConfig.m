function params = preConfig(config)

config = lower(config);

switch config
    
    case 'default'
        
        %DEFAULT GENERAL PARAMETERS
        params.preConfig = 'default';
        params.fs = 1;                          %numeric [default=1]
        params.resetSeed = true;                %true[default] or false
        params.openParPool = true;              %true[default] or false
        
        %DEFAULT EMD PARAMETERS
        params.EMDStopThresh = 5;               %numeric [default=5]
        params.keepDC = false;                  %true or false[default]
        params.computeNormalized = true;        %true[default] or false
        params.EthreshPCT = 99.999;             %numeric [default=99.999]
        params.VerboseEMD = false;              %true or false[default]
        params.statusEMD = true;                %true[default] or false
        
        %DEFAULT CEEMD PARAMETERS
        params.trials = 5;                          %numeric [default=??????????]
        params.maskSignal = 'tone';                 %'tone'[default] or 'filteredNoise' or 'siftedNoise'
        params.maskSignalcutoff = 'IAweightedIF';   %'IAweightedIF'[default] or 'EweightedIF' or 'meanIF'
        params.maskAmplitude = 3;                   %numeric [default=3]
        params.maskFreqSpacing = 0.75;              %numeric [default=0.75]
        params.maskNoiseType = 'normal';            %'normal'[default] or 'uniform'
        params.maskFilterOrder = 12;                %numeric [default=12]
        
        %DEFAULT DEMODLATION PARAMETERS
        params.demodMethod = 'VMDdemod';        %'VMDdemod'[default] or 'IMFdemod' or 'HTdemod'
        params.IAestMaxIter = 5;                %numeric [default=5]
        params.IFestThrowOut = 0;               %numeric [default=0]
        params.VMDdemodAlpha = 5e4;             %numeric [default=5e4]
        params.VMDdemodTau = 1e-4;              %numeric [default=1e-4]
        params.IMFsigmaEstDeriv = 'center3';    %see derivApprox
        params.demodPhaseDeriv = 'center9';     %see derivApprox
        params.demodFiltFreq = 1/150;           %numeric or Inf (for none)
        
        %DEFAULT SIFTING PARAMETERS
        params.SiftStopThresh = -120;           %numeric [default=-120]
        params.maxSifts = 500;                  %numeric [default=500]
        params.siftAlpha = 1;                   %numeric [default=1]
        params.VerboseSift = false;             %true or false
        params.ExtremaID = 'parabolic';         %'parabolic' or 'simple'
        params.postextrapolate = 'oddReflect';  %'oddReflect' or 'none'
        params.tangentInterpIter = 0;           %numeric [default=0]
        params.extremaEndpoints = 'discard';    %'discard'[default] or 'keep'
        params.lpcExtrapLen = 200;              %numeric [default=200]
        params.numLPC = 200;                    %numeric [default=200]
        params.tangentDeriv = 'center9';        %see derivApprox [default='center9']
        params.siftInterp = 'spline';           %'spline'[default]
        params.zeroMeanIMFs = true;             %true[default] or false
        
    case 'emd'
        
        %DEFAULT GENERAL PARAMETERS
        params.preConfig = 'emd';
        params.fs = 1;                          %numeric [default=1]
        params.resetSeed = true;                %true[default] or false
        params.openParPool = false;             %true[default] or false
        
        %DEFAULT EMD PARAMETERS
        params.EMDStopThresh = 5;               %numeric [default=5]
        params.keepDC = false;                  %true or false[default]
        params.computeNormalized = true;        %true[default] or false
        params.EthreshPCT = 100;                %numeric [default=99.999]
        params.VerboseEMD = false;              %true or false[default]
        params.statusEMD = true;                %true[default] or false
        
        %DEFAULT CEEMD PARAMETERS
        params.trials = 1;                          %numeric [default=??????????]
        params.maskSignal = 'tone';                 %'tone'[default] or 'filteredNoise' or 'siftedNoise'
        params.maskSignalcutoff = 'IAweightedIF';   %'IAweightedIF'[default] or 'EweightedIF' or 'meanIF'
        params.maskAmplitude = 0;                   %numeric [default=3]
        params.maskFreqSpacing = 0.75;              %numeric [default=0.75]
        params.maskNoiseType = 'normal';            %'normal'[default] or 'uniform'
        params.maskFilterOrder = 12;                %numeric [default=12]
        
        %DEFAULT DEMODLATION PARAMETERS
        params.demodMethod = 'VMDdemod';        %'VMDdemod'[default] or 'IMFdemod' or 'HTdemod'
        params.IAestMaxIter = 5;                %numeric [default=5]
        params.IFestThrowOut = 0;               %numeric [default=0]
        params.VMDdemodAlpha = 5e4;             %numeric [default=5e4]
        params.VMDdemodTau = 1e-4;              %numeric [default=1e-4]
        params.IMFsigmaEstDeriv = 'center3';    %see derivApprox
        params.demodPhaseDeriv = 'center9';     %see derivApprox
        params.demodFiltFreq = 1/150;           %numeric or Inf (for none)
        
        %DEFAULT SIFTING PARAMETERS
        params.SiftStopThresh = -120;           %numeric [default=-120]
        params.maxSifts = 500;                  %numeric [default=500]
        params.siftAlpha = 1;                   %numeric [default=1]
        params.VerboseSift = false;             %true or false
        params.ExtremaID = 'simple';            %'parabolic' or 'simple'
        params.postextrapolate = 'none';        %'oddReflect' or 'none'
        params.tangentInterpIter = 0;           %numeric [default=0]
        params.extremaEndpoints = 'keep';       %'discard'[default] or 'keep'
        params.lpcExtrapLen = 0;                %numeric [default=200]
        params.numLPC = 200;                    %numeric [default=200]
        params.tangentDeriv = 'center9';        %see derivApprox [default='center9']
        params.siftInterp = 'spline';           %'spline'[default]
        params.zeroMeanIMFs = false;             %true[default] or false
        
    case 'hht'
        
        params.preConfig = 'hht';
        params.fs = 1;                          %numeric [default=1]
        params.resetSeed = true;                %true[default] or false
        params.openParPool = false;             %true[default] or false
        
        %DEFAULT EMD PARAMETERS
        params.EMDStopThresh = 5;               %numeric [default=5]
        params.keepDC = false;                  %true or false[default]
        params.computeNormalized = false;       %true[default] or false
        params.EthreshPCT = 100;                %numeric [default=99.999]
        params.VerboseEMD = false;              %true or false[default]
        params.statusEMD = true;                %true[default] or false
        
        %DEFAULT CEEMD PARAMETERS
        params.trials = 1;                          %numeric [default=??????????]
        params.maskSignal = 'tone';                 %'tone'[default] or 'filteredNoise' or 'siftedNoise'
        params.maskSignalcutoff = 'IAweightedIF';   %'IAweightedIF'[default] or 'EweightedIF' or 'meanIF'
        params.maskAmplitude = 0;                   %numeric [default=3]
        params.maskFreqSpacing = 0.75;              %numeric [default=0.75]
        params.maskNoiseType = 'normal';            %'normal'[default] or 'uniform'
        params.maskFilterOrder = 12;                %numeric [default=12]
        
        %DEFAULT DEMODLATION PARAMETERS
        params.demodMethod = 'HTdemod';         %'VMDdemod'[default] or 'IMFdemod' or 'HTdemod'
        params.IAestMaxIter = 5;                %numeric [default=5]
        params.IFestThrowOut = 0;               %numeric [default=0]
        params.VMDdemodAlpha = 5e4;             %numeric [default=5e4]
        params.VMDdemodTau = 1e-4;              %numeric [default=1e-4]
        params.IMFsigmaEstDeriv = 'forward';    %see derivApprox
        params.demodPhaseDeriv = 'forward';     %see derivApprox
        params.demodFiltFreq = Inf;             %numeric or Inf (for none)
        
        %DEFAULT SIFTING PARAMETERS
        params.SiftStopThresh = -120;           %numeric [default=-120]
        params.maxSifts = 500;                  %numeric [default=500]
        params.siftAlpha = 1;                   %numeric [default=1]
        params.VerboseSift = false;             %true or false
        params.ExtremaID = 'simple';            %'parabolic' or 'simple'
        params.postextrapolate = 'none';        %'oddReflect' or 'none'
        params.tangentInterpIter = 0;           %numeric [default=0]
        params.extremaEndpoints = 'keep';       %'discard'[default] or 'keep'
        params.lpcExtrapLen = 0;                %numeric [default=200]
        params.numLPC = 200;                      %numeric [default=200]
        params.tangentDeriv = 'center9';        %see derivApprox [default='center9']
        params.siftInterp = 'spline';           %'spline'[default]
        params.zeroMeanIMFs = false;             %true[default] or false
        
    case 'ceemdanht'
        
        %DEFAULT GENERAL PARAMETERS
        params.preConfig = 'ceemdan';
        params.fs = 1;                          %numeric [default=1]
        params.resetSeed = true;                %true[default] or false
        params.openParPool = true;              %true[default] or false
        
        %DEFAULT EMD PARAMETERS
        params.EMDStopThresh = 5;               %numeric [default=5]
        params.keepDC = false;                  %true or false[default]
        params.computeNormalized = true;        %true[default] or false
        params.EthreshPCT = 99.999;             %numeric [default=99.999]
        params.VerboseEMD = false;              %true or false[default]
        params.statusEMD = true;                %true[default] or false
        
        %DEFAULT CEEMD PARAMETERS
        params.trials = 50;                         %numeric [default=50]
        params.maskSignal = 'siftedNoise';          %'tone'[default] or 'filteredNoise' or 'siftedNoise'
        params.maskSignalcutoff = 'IAweightedIF';   %'IAweightedIF'[default] or 'EweightedIF' or 'meanIF'
        params.maskAmplitude = 3;                   %numeric [default=3]
        params.maskFreqSpacing = 0.75;              %numeric [default=0.75]
        params.maskNoiseType = 'normal';            %'normal'[default] or 'uniform'
        params.maskFilterOrder = 12;                %numeric [default=12]
        
        %DEFAULT DEMODLATION PARAMETERS
        params.demodMethod = 'HTdemod';         %'VMDdemod'[default] or 'IMFdemod' or 'HTdemod'
        params.IAestMaxIter = 5;                %numeric [default=5]
        params.IFestThrowOut = 0;               %numeric [default=0]
        params.VMDdemodAlpha = 5e4;             %numeric [default=5e4]
        params.VMDdemodTau = 1e-4;              %numeric [default=1e-4]
        params.IMFsigmaEstDeriv = 'forward';    %see derivApprox
        params.demodPhaseDeriv = 'forward';     %see derivApprox
        params.demodFiltFreq = Inf;           %numeric or Inf (for none)
        
        %DEFAULT SIFTING PARAMETERS
        params.SiftStopThresh = -120;           %numeric [default=-120]
        params.maxSifts = 500;                  %numeric [default=500]
        params.siftAlpha = 1;                   %numeric [default=1]
        params.VerboseSift = false;             %true or false
        params.ExtremaID = 'simple';            %'parabolic' or 'simple'
        params.postextrapolate = 'none';  %'oddReflect' or 'none'
        params.tangentInterpIter = 0;           %numeric [default=0]
        params.extremaEndpoints = 'keep';    %'discard'[default] or 'keep'
        params.lpcExtrapLen = 0;              %numeric [default=200]
        params.numLPC = 200;                    %numeric [default=200]
        params.tangentDeriv = 'center9';        %see derivApprox [default='center9']
        params.siftInterp = 'spline';           %'spline'[default]
        params.zeroMeanIMFs = false;             %true[default] or false
        
        
    otherwise
        error('preConfig: invalid configuration type.')
end