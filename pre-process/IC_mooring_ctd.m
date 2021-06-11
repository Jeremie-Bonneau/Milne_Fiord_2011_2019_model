function [SIC_mesh,TIC_mesh] = IC_mooring_ctd(empty_mesh,t0,y_top,y_bot,dy)
% SIC_mooring writes the salinity initial conditions by fitting an tanh
% profile with the salinity data available from the mooring and also the
% temperature IC with the mooring data.
%   
%INPUTS
%    empty_mesh : nan mesh created by meshing
%    t0         : Initial time                                     [serial]
%    y_top      : Top depth                                             [m]
%    y_bot      : Bottom depth                                          [m]
%    dy         : Vertical discretization                               [m]

%OUTPUT
%    SIC_mesh   : mesh with the absolute salinity IC on the first column
%    TIC_mesh   : mesh with the consevative temperature IC on the first
%                 column

%VERSION 1, August 2019

%% DISPLAY
disp('IC_mooring_ctd : Writing initial conditions')

%% VERIFICATIONS

% Check number of input
if ~(nargin == 5)   %Checking number of args
    error('IC_mooring_ctd : Requires 5 inputs')
end %if

% Check time input
if ~(isnumeric(t0))
    error('IC_mooring_ctd : t0 has to be in serial time')
end

% Check mesh vertical discretization
l=length(empty_mesh(:,1));
if (y_bot-y_top)/dy~=(l-1)
    error('IC_mooring_ctd : input mesh not the right size')
end

%% TEMPERATURE DATA
% [T,TY,Time] = MEL_Tmatrix('TBD',t0,t0+1,'day',0,40,0.1,'linear');
% if t0-Time(1)>2
%     disp('IC_mooring_ctd : WARNING! , Mooring data > 2 days away from t0');
% end
% T=T(:,1);
Y=[y_top:dy:y_bot]';
ds=importdata('data\Tmatrix_0_30_day_10cm.mat'); % Import 30min avg data
T=interp2(ds.Time,ds.Y,ds.Tmatrix,t0,Y);

%% SALINITY DATA
ds=importdata('data\SA_IC_ctd.mat');
%%
y=year(datetime(t0,'ConvertFrom','datenum'));
if y==2011
    ds=ds(1);
elseif y==2012
    ds=ds(2);
elseif y==2013
    ds=ds(3);
elseif y==2014
    ds=ds(4);
elseif y==2015
    ds=ds(5);
elseif y==2016
    ds=ds(6);
elseif y==2017
    ds=ds(7);
elseif y==2018
    ds=ds(8);
else
    error('IC_mooring_ctd : t0 is not inbound')
end

if ds.date~=t0
    disp('IC_mooring_ctd : WARNING! t0 is not the right day for ctd IC')
end

%% WRITING IC

% Salt
SIC_mesh=empty_mesh;
SIC_mesh(:,1)=interp1(ds.press,ds.SA,Y);

% Temperature
TIC_mesh=empty_mesh;
TIC_mesh(:,1)=gsw_CT_from_t(SIC_mesh(:,1),T,Y);

end