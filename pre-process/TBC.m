function [TBC_mesh] = TBC(TIC_mesh,t0,tf,dt,y_top,y_bot)
%TBC writes the temperature boundary conditions on the mesh based on
%the mooring 30min averaged matrix data (Tmatrix_0_30_30m_10cm). Averaging
%when dt > 60 minutes.
%   
%INPUTS
%    -  TIC_mesh  : nan mesh with initial conditions on the first column
%    -  t0        : Initial time                                   [serial]
%    -  tf        : Final time                                     [serial]
%    -  dt        : Time discretization                               [min]
%    -  y_top     : Top depth                                           [m]
%    -  y_bot     : Bottom depth                                        [m]
%
%OUTPUT
%    -  TBC_mesh  : Mesh with the boundary conditions on the first and last
%                   rows (and inital conditions on first column)
% 
%VERSION 1.1, August 2019

%% DISPLAY
disp('TBC : writing temperature boundary conditions')

%% VERIFICATIONS

% Check number of input
if ~(nargin == 6)   %Checking number of args
    error('TBC :  Requires 6 inputs')
end %if

% Check y_top and Y_bot
if ~mod(y_top,0.1)==0
    error('TBC : top of mesh has to be a 10cm multiple')
end
if ~mod(y_bot,0.1)==0
    error('TBC : top of mesh has to be a 10cm multiple')
end

% Check time
if dt>=60
    if ~mod(dt,30)==0
        error('TBC : time interval has to be a multiple of 30 min')
    end
end

%% TEMPERATURE DATA

ds=importdata('data\Tmatrix_0_30_30m_10cm.mat'); % Import 30min avg data
%% Top
[~,i_top]=min(abs(ds.Y-y_top));     %find top index
% Bottom
[~,i_bot]=min(abs(ds.Y-y_bot));     %find bottom index
% Time
[~,i_0]=min(abs(ds.Time-t0));    %find initial time index
[~,i_f]=min(abs(ds.Time-tf));    %find final time index

% Resampling or downsampling
time=t0:1/(24*60/dt):tf; %model matrix time array
if dt<60
    TIC_mesh(1,:)=interp1(ds.Time(i_0-1:i_f+1),ds.Tmatrix(i_top,i_0-1:i_f+1),time);   %resampling
    TIC_mesh(end,:)=interp1(ds.Time(i_0-1:i_f+1),ds.Tmatrix(i_bot,i_0-1:i_f+1),time); %resampling  
else
    n=dt/30; %number of samples to average
    T_top=ds.Tmatrix(i_top,i_0-n:i_f+2*n);      %slicing top temp
    T_bot=ds.Tmatrix(i_bot,i_0-n:i_f+2*n);      %slicing bottom temp
    T_time=ds.Time(i_0-n:i_f+2*n);              %slicing time
    M=length(T_top)-mod(length(T_top),n);

    T_top=sum(reshape(T_top(1:M),n,[]))./n;     %downsampling top
    T_bot=sum(reshape(T_bot(1:M),n,[]))./n;     %downsampling bottom
    T_time=sum(reshape(T_time(1:M),n,[]))./n;   %downsampling time
    
    TIC_mesh(1,:)=interp1(T_time,T_top,time);   %resampling top
    TIC_mesh(end,:)=interp1(T_time,T_bot,time); %resampling bottom
end

%% Conservative temperature
% first pass approx
if y_top<10
    SA_top=zeros(1,length(TIC_mesh(1,:)))+0.1;
else
    SA_top=zeros(1,length(TIC_mesh(1,:)))+27.0;
end
if y_bot<10
    SA_bot=zeros(1,length(TIC_mesh(1,:)))+0.1;
else
    SA_bot=zeros(1,length(TIC_mesh(1,:)))+27.0;
end

TIC_mesh(1,:)=gsw_CT_from_t(SA_top,TIC_mesh(1,:),y_top);
TIC_mesh(end,:)=gsw_CT_from_t(SA_bot,TIC_mesh(end,:),y_bot);

TBC_mesh=TIC_mesh;
TIC_mesh=[];
    
end