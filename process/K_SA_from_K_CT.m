function [K_SA] = K_SA_from_K_CT(K_CT)
% K_SA_from_K_CT computes the diffusivity (molecular+turb) for salinity
% from the turbulent eddy coefficiwent for heat using the parameterization
% proposed by Jackson and Rehman (2014)

% The parameter eps/nu/N^2 is first estimated with the formula from table 2
% using the SC method (see Jackson and Rehman 2014):
%
%              K_CT/nu = a(y/y_t)^n1     if y<y_t 
%                      = a(y/y_t)^n2     if y>y_t
% y=eps/nu/N^2, y_t=0.36, a=0.31, n1=0.16, n2=1.06
%
% Then the ratio K_SA/K_CT is computed using equation 6:
%
%        1+K_m_SA/K_m_CT     1-K_m_SA/K_m_CT
%    d = _______________ + ( _______________ ) * tanh(a1[log10(y)-a2])
%               2                   2
% a1=0.92, a2=0.60, K_m_SA/K_m_CT is 0.01

%INPUTS
%    K_CT     : Heat diffusivity (mol+turb)                      [m^2 s^-1]
%
%OUTPUT
%    K_SA     : Salt diffusivity (mol+turb)                      [m^2 s^-1]
%
%VERSION 1, Sept 2019

%% DISPLAY
%disp('K_SA_from_K_CT : Computing K_SA')

%% VERIFICATIONS

% Check number of input
if ~(nargin == 1)   %Checking number of args
    error('K_SA_from_K_CT :  Requires 1 input')
end 

%% loading structure;
ds=importdata('data\K_SA_from_K_CT_struct2interp.mat');

%% interpolation from structure arrays
K_SA=interp1(ds.K_CT,ds.K_SA,K_CT);

end
