function [A,IF,PSI,SIGMA] = HTdemod(S,varargin)

params.demodPhaseDeriv = 'center9';
params.fs = 1;

params = parse_pv_pairs(params,varargin);




PSI = hilbert(S);
[A,IF,~,SIGMA] = amfmdemod(PSI,'fs',params.fs,'demodPhaseDeriv',params.demodPhaseDeriv);