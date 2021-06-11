function [CT_diff,CT_mooring,stats] = CT_Eval_fast(t0,tf,dt,y_top,y_bot,dy,SAmatrix,CTmatrix,TmatrixEval)
% CT_Eval_fast helps evaluating the model by computing the conservative
% temperature difference between the model and the thermistors data. 
% Because of the irregular thermistor's spacing, linear interpolation is 
% used between the thermistors so it kinda weights each thermistors data.
% Absolute salinity returned by the model is used to compute conservative
% temperature from the mooring data. fast version uses the evaluation
% temperature matrix as an input instead for building it from the database
%   
%   
%INPUTS
%    -  t0        : Initial time                                   [serial]
%    -  tf        : Final time                                     [serial]
%    -  dt        : Time discretization                               [min]
%    -  y_top     : Top depth                                           [m]
%    -  y_bot     : Bottom depth                                        [m]
%    -  dy        : Vertical discretization                             [m]
%    -  CTmatrix  : Temperature matrix returned by the model         [degC]
%    -  SAmatrix  : Absolute salinity returned by the model          [g/kg]
%    -  TmatrixEval: Temperature matrix structure used for evauation
%                   contains Tmatrix (mxn), Y (1xm), Time (1xn).
%
%OUTPUT
%    -  CT_diff   : Conservative temperature difference matrix, positive
%                   if the model overevaluates                       [degC]
%    -  CT_mooring: Conservative temperature mooring data            [degC]
%    -  stats     : Table with stats for the whole time

%VERSION 1, August 2019

%% DISPLAY
%%disp('Conservative Temperature evaluation')

%% Verifications
% Check number of input
if ~(nargin == 9)   %Checking number of args
    error('CT_Eval_fast :  Requires 9 inputs')
end

% Check time inputs
if ~(isnumeric(t0)&&isnumeric(tf))
    error('CT_Eval_fast : t0 and tf have to be in serial time')
end
if mod(((tf-t0)*24*60),dt)~=0
    error('CT_Eval_fast : (tf-t0)/dt has a remainder')
end

% Check spacial input
if mod((y_bot-y_top),dy)~=0
    error('CT_Eval_fast : (y_bot-y_top)/dt has a remainder')
end

%% Check Tmatrix dimensions
[m,n]=size(CTmatrix);
if m~=((y_bot-y_top)/dy)+1
    error('CT_Eval_fast : Tmatrix not the right size')
end
if n~=((tf-t0)*24*60/dt)+1
    error('CT_Eval_fast : Tmatrix not the right size')
end
if size(CTmatrix)~=size(SAmatrix)
    error('CT_Eval_fast : CTmatrix and SAmatrix not the same size');
end

%% Thermistor data
ds=TmatrixEval; % Import 30 min avg data
%%
th_matrix=interp1(ds.Time',ds.Tmatrix',[t0:dt/60/24:tf]');
th_matrix=interp1(ds.Y',th_matrix',[y_top:dy:y_bot]');
P=repmat([y_top:dy:y_bot]',1,n);
th_CT=gsw_CT_from_t(SAmatrix,th_matrix,P);
CT_mooring=th_CT;

%% Evaluation
CT_diff=CTmatrix-th_CT;

Mean=mean(CT_diff,'all');
Median=median(CT_diff,'all');
STD=std(CT_diff,0,'all');
Max=max(CT_diff,[],'all');
Min=min(CT_diff,[],'all');
RMSE=rms(reshape(CT_diff,1,m*n));
MAE=mean(abs(CT_diff),'all');
SStot=sum((th_CT-mean(th_CT,'all')).^2,'all');
SSres=sum((CT_diff).^2,'all');
R2=1-SSres/SStot;
Values=[Mean; Median; STD; Max; Min; RMSE; MAE; R2];

stats=table(Values);
stats.Properties.RowNames={'Mean';'Median';'STD';'Max';'Min';'RMSE';'MAE';'R2'};

end