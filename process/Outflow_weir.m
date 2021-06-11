function [SA_profile_out,CT_profile_out,outflow] = Outflow_weir(SA_profile_in,CT_profile_in,N2,y_top,y_bot,dy,coef,h0)
% Outflow simulates the outflow of water from the epishelf lake. It pulls
% out the specified amount of water between h0 and the base of the
% halocline and resamples the water profile adding the same amount of water
% at the base of the profile using the same water properties than the
% Profile_in base node. 
%
% It uses the Kindsvater-Carter weir equation for rectangular weir
%
%                   Q = 2/3 * (2g)^0.5 * C * b * h^1.5           [m^3 s^-2]
%
% where C is coef for friction, b is the width of the weir and h is the
% height of water above the weir. Replacing the gravity by the reduced
% gravity (g') of the system (g*delta rho/rho ~ 0.245), C*b by coef and 
% defining h as the depth of outflowing water between the minimum ice shelf
% draft (h0) and the bottom of halocline (N2=0.01) we get:
% 
%                dh/dt =  2/3 * (2g')^0.5 * Coef * h^1.5 / A       [m s^-1]
% 
% Where A is the area of the lake (~72 km^2 = 72x10^6 m^2)

%   
%INPUTS
%    SA_Profile_in : Conservative temperature profile before outflow [degC]
%    CT_profile_in : Absolute salinity profile before outlfow        [g/kg]
%    N2            : N2 profile                                      [s^-2]
%    y_top         : top of the water column under sea surface          [m]
%    y_bot         : bottom of the water column under sea surface       [m]
%    dy            : Vertical discretization                            [m] 
%    coef          : friction coefficient * width of channel            [m]
%    h0            : top depth of the outflow                           [m]
%
%OUTPUT
%    SA_Profile_out: Conservative temperature profile after outflow  [degC]
%    CT_profile_out: Absolute salinity profile after outlfow         [g/kg]
%
%VERSION 1, sept 2019
%
%% DISPLAY
fprintf('Outflow_weir , h0=%2.2f , coef=%2.2f \n',h0,coef)

%% VERIFICATIONS

% Check number of input
if ~(nargin == 8)   %Checking number of args
    error('Outflow_weir :  Requires 8 inputs')
end

% Check profiles
[~,n1]=size(SA_profile_in);
[~,n2]=size(CT_profile_in);

if sum(isnan(CT_profile_in),'all')~=0
    error('Outflow_weir : error in CT_profile: has NaN value(s)')
end
if sum(isnan(SA_profile_in),'all')~=0
    error('Outflow_weir : error in SA_profile: has NaN value(s)')
end
if size(CT_profile_in)~=size(SA_profile_in)
    error('Outflow_weir : error in profiles sizes, do not match')
end
if sum(isnan(N2),'all')~=0
    error('Outflow_weir : error in N2 profile: has NaN value(s)')
end
if size(N2)~=size(CT_profile_in)
    error('Outflow_weir : N2 profile size do not match')
end
if n1~=1
    error('Outflow_weir : SA_profile_in has to be a column array')
end
if n2~=1
    error('Outflow_weir : CT_profile_in has to be a column array')
end

% Check spacial input
if mod((y_bot-y_top),dy)~=0
    error('Outflow_weir : (y_bot-y_top)/dy has a remainder')
end
if mod(dy,0.001)~=0
    error('Outflow_weir : dy has to be a 0.001 multiple')
end

%% OUTFLOW DEPTH
Y=[y_top:dy:y_bot]';            %model vertical array
top=h0;                         %depth at top of outflow layer
[~,idx_top]=min((Y-top).^2);    %index top outflow depth
idx=find(N2>0.01);
idx_bot=idx(end);               %index of bottom of halocline
bot=Y(idx_bot);                 %depth of bottom of halocline
h=bot-top;                      %outflow layer thickness
m=(idx_bot-idx_top);            %number of dy in the outflow layer

%% CALCULATE DAY OUTFLOW
outflow = 2/3 * (2*0.245)^0.5 * coef * h^1.5 / 72000000 * 24*60*60;
outflow=round(outflow,3);
n=outflow/0.001;
n=round(n,0);
%% RESAMPLING PROFILE at 1mm
r=dy/0.001;
Y2=[y_top:0.001:y_bot]';
[~,idx_top]=min((Y2-top).^2);  %top outflow layer index
[~,idx_bot]=min((Y2-bot).^2);  %bottom outflow layer index
SA_profile2=interp1(Y,SA_profile_in,Y2);
CT_profile2=interp1(Y,CT_profile_in,Y2);
Yout=(linspace(Y2(idx_top),Y2(idx_bot),m*r+1-n))'; %shrinked halocline depths

%% SA
SAtop=SA_profile2(1:idx_top-1);
SAout=SA_profile2(idx_top:idx_bot);
SAbot=SA_profile2(idx_bot+1:end);
SAout=interp1((Y2(idx_top):0.001:Y2(idx_bot))',SAout,Yout); %shrinked halocline
SA_profile2=[SAtop;SAout;SAbot];

%% CT
CTtop=CT_profile2(1:idx_top-1);
CTout=CT_profile2(idx_top:idx_bot);
CTbot=CT_profile2(idx_bot+1:end);
CTout=interp1((Y2(idx_top):0.001:Y2(idx_bot))',CTout,Yout); %shrinked halocline
CT_profile2=[CTtop;CTout;CTbot];

%% BASAL INFLOW
SAin=repmat(SA_profile2(end),n,1);
CTin=repmat(CT_profile2(end),n,1);
SA_profile2=[SA_profile2;SAin];     %adding same water amount at the base of the profile
CT_profile2=[CT_profile2;CTin];     %adding same water amount at the base of the profile

%% RESAMPLING AT dy
SA_profile_out=interp1(Y2,SA_profile2,Y);
CT_profile_out=interp1(Y2,CT_profile2,Y);


end