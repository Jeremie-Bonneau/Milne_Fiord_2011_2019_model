function [mesh] = meshing(t0,tf,dt,y_top,y_bot,dy)
% meshing creates a mesh (space in columns and time in rows)

%INPUT
%   t0 is the starting time, serial
%   tf is the final time, serial
%   dt is the time step, in minutes
%   y_top is the depth at the top, in meter
%   Y_bot is the depth at the bottom, in meter
%   dy is the vertical discretization, in meter

%OUTPUT
%   mesh: nan matrix with appropriate dimensions

%VERSION 1, May 2019

%% DISPLAY
disp('meshing : Creating mesh')

%% VERIFICATIONS

% Check number of input
if ~(nargin == 6)   %Checking number of args
    error('mesh:  Requires 6 inputs')
end %if

% Check time inputs
if ~(isnumeric(t0)&&isnumeric(tf))
    error('mesh: t0 and tf have to be in serial time')
end
if mod(((tf-t0)*24*60),dt)~=0
    error('mesh: (tf-t0)/dt has a remainder')
end

% Check spacial input
if mod((y_bot-y_top),dy)~=0
    error('mesh: (y_bot-y_top)/dt has a remainder')
end

%% CREATE MESH

m=(y_bot-y_top)/dy+1;   %number of rows
n=((tf-t0)*24*60/dt)+1; %number of column
mesh=nan(m,n);

end