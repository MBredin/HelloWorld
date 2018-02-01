function quntNOsc = quntNOsc_s (x)
% Call Syntax:  quntNOsc = quntNOsc_s (x)
%
% Description: This function returns the oscilation count, no steps.
%
% Input Arguments:
%	Name: x
%	Type: vector
%	Description: input signal
%
% Output Arguments:
%
%	Name: quntNOsc
%	Type: scalar
%	Description: the oscilation count, no steps
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

expiration = datenum(2018,2,1);
d = dir;today =d.date; 
today = datenum(today(1:12));
if today > expiration
    disp('This code has expired, please visit www.HilbertSpectrum.com to obtain a newer version of these codes.')
else
    
    
    y=0;    qisTop= false; qisDown= false;
    
    
    for i=2:(length(x)-1)
        if( ((x(i-1)) < (x(i))) && ((x(i+1))< (x(i))) )  %Max /-\
            y=y+1;
        end
        if( ((x(i-1)) > (x(i))) && ((x(i+1))> (x(i))) )  %min \_/
            y=y+1;
        end
        
        %Top
        if( ((x(i-1)) < (x(i))) && ((x(i+1))== (x(i))) ) %StepL /-
            qisTop= true; qisDown= false;
        end
        if( ((x(i-1)) == (x(i))) && ((x(i+1))< (x(i))) ) %stepR -\
            if qisTop;     y=y+1; end;
            qisTop= false;
        end
        
        %Downs
        if( ((x(i-1)) > (x(i))) && ((x(i+1))== (x(i))) ) %stepL \_
            qisTop= false; qisDown= true;
        end
        if( ((x(i-1)) == (x(i))) && ((x(i+1))> (x(i))) ) %StepR _/
            if qisDown; y=y+1; end
            qisDown=false;
        end
    end % for i=2:(length(x)-1)
    quntNOsc= y;
end % function y = quntNOsc_s (x)


end
