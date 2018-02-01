function [u_pAndt_p, l_qAndt_q] = getExtrama(x,varargin)

%OPTIONS
%   method = 'parabolic' [default] or 'simple'
%   extrapolate = 'oddReflect' [default] or 'none'

%------------------
% Check valid input
%------------------

params.method = 'parabolic';
params.postextrapolate = 'oddReflect';
params.extremaEndpoints = 'discard';
params = parse_pv_pairs(params,varargin);

%-----------
% Initialize
%-----------


%-----
% Main
%-----

%FIND EXTREMA
if strcmp(params.method, 'parabolic')  %PARABOLIC EXTREMA
    Maxs = parabolicMaxs(x,params.extremaEndpoints);    %get parabolic maxima
    Mins = parabolicMins(x,params.extremaEndpoints);    %get parabolic minima
else                            %SIMPLE EXTREMA
    Maxs = simpleMaxs(x,params.extremaEndpoints);       %get simple maxima
    Mins = simpleMins(x,params.extremaEndpoints);       %get simple minima
end

%EXTRAPOLATE EXTREMA 
if strcmp(params.postextrapolate, 'oddReflect')                %EXTRAPOLATION                            
    u_pAndt_p = extrapolMaxs(Maxs, Mins, length(x));    %extrapolate maxima (prevent end effects)
    l_qAndt_q = extrapolMins(Maxs, Mins, length(x));    %extrapolate minima (prevent end effects)
elseif strcmp(params.postextrapolate, 'none')                  %NO EXTRAPOLTION
    u_pAndt_p = Maxs;                                   %maxima
    l_qAndt_q = Mins;                                   %minima
else
    error('getExtrama: unknown extrema extrapolation method')
end






