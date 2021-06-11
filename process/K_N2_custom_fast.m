function [Ks,Kh] = K_N2_custom_fast(N2,y_top,y_bot,dy,K_CT_top,K_CT_bot)
% K_N2_custom computes the combined (molecular+turbulent) mixing for heat
% and salt using the N2 frequency profile. 
%
%INPUTS
%    N2         : N2 profile                                         [s^-2]
%    y_top      : top of the water column under sea surface             [m]
%    y_bot      : bottom of the water column under sea surface          [m]
%    dy         : Vertical discretization                               [m] 
%    K_CT_top   : Heat diffusivity coef (mol+turb) above halocline  [m^2/s]
%    K_CT_bot   : Heat diffusivity coef (mol+turb) below halocline  [m^2/s]

%OUTPUT
%    Kh         : Heat diffusivity coef (mol+turb) column           [m^2/s]
%    Ks         : Salt diffusivity coef (mol+turb) column           [m^2/s]
%
%VERSION 2, Sept 2019

%% DISPLAY
%disp('K_N2_custom_fast : Computing mixing coefficients')

%% VERIFICATIONS
% 
% % Check number of input
% if ~(nargin == 6)   %Checking number of args
%     error('K_N2_custom :  Requires 6 inputs')
% end %if
% 
% % Check spacial input
% if mod((y_bot-y_top),dy)~=0
%     error('K_N2_custom : (y_bot-y_top)/dt has a remainder')
% end
% 
% % Check N2 profile
% if sum(isnan(N2),'all')~=0
%     error('K_N2_custom : error N2 profile: has NaN value(s)')
% end
% if size(N2)~=size([y_top:dy:y_bot]')
%     error('K_N2_custom : error in profiles sizes, do not match')
% end

%% SOME VARIABLES
Y=[y_top:dy:y_bot]';                 %vertical array
[~,N2i]=max(N2);                     %max N2 index
Ymax=Y(N2i);                         %depth at which N2 is max

%% FIND K_SA
K_SA_top=K_SA_from_K_CT(K_CT_top);
K_SA_bot=K_SA_from_K_CT(K_CT_bot);


%% Heat diffusion
Kh=zeros(size(Y));                   %creating array
K_CT_m=1.4*10^(-7);                  %thermal molecular diffusivity [m^2/s]
Kh=Kh+[N2>=0.01]*K_CT_m;             %add molecular diffusion at halocline
Kh=Kh+[N2<0.01 & Y<Ymax]*K_CT_top;   %top layer diffusivity
Kh=Kh+[N2<0.01 & Y>Ymax]*K_CT_bot;   %bottom layer diffusivity

%% Salt diffusion
Ks=zeros(size(Y));                   %creating array
K_SA_m=1.4*10^(-9);                  %salt molecular diffusivity    [m^2/s]
Ks=Ks+[N2>=0.01]*K_SA_m;             %add molecular diffusion at halocline
Ks=Ks+[N2<0.01 & Y<Ymax]*K_SA_top;   %top layer diffusivity
Ks=Ks+[N2<0.01 & Y>Ymax]*K_SA_bot;   %bottom layer diffusivity

end