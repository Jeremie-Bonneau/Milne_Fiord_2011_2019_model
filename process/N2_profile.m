function [N2] = N2_profile(SA_profile,CT_profile,y_top,y_bot,dy)
% N2 computes the N2 frequency of the water column using the gsw toolbox 
%   
%INPUTS
%    CT_profile : Conservative temperature profile                   [degC]
%    SA_profile : Absolute salinity profile                          [g/kg]
%    y_top      : top of the water column under sea surface             [m]
%    y_bot      : bottom of the water column under sea surface          [m]
%    dy         : Vertical discretization                               [m] 
%
%OUTPUT
%    N2         : N2 profile                                         [s^-2]
%
%VERSION 1, August 2019
%
%% DISPLAY
disp('N2 : Computing N2')

%% VERIFICATIONS

% Check number of input
if ~(nargin == 5)   %Checking number of args
    error('N2 :  Requires 5 inputs')
end

% Check profiles
if sum(isnan(CT_profile),'all')~=0
    error('N2 : error in CT_profile: has NaN value(s)')
end
if sum(isnan(SA_profile),'all')~=0
    error('N2 : error in SA_profile: has NaN value(s)')
end
if size(CT_profile)~=size(CT_profile)
    error('N2 : error in profiles sizes, do not match')
end

% Check spacial input
if mod((y_bot-y_top),dy)~=0
    error('N2 : (y_bot-y_top)/dy has a remainder')
end

%% N2 COMPUTATION
[N2,y_mid]=gsw_Nsquared(SA_profile,CT_profile,[y_top:dy:y_bot]',80);
N2=interp1(y_mid,N2,[y_top:dy:y_bot]');   % interpolation at model nodes
N2(1)=N2(2);   N2(end)=N2(end-1);         % patching first and last values nan

end