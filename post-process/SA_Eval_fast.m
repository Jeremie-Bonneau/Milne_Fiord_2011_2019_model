function [SA_diff,SA_mooring,Y,stats] = SA_Eval_fast(t0,tf,dt,y_top,y_bot,dy,SAmatrix,SA_eval)
% SA_Eval helps evaluating the model by computing the absolute salinity 
% difference between the model and the conductivity measurements on the 
% mooring. The salinity at each instrument is compared with the model
% output. Only works for a single mooring deployment.
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
%    -  SA_eval   : Struture with 3 fields, SAmatrix (mxn), Y (mxn) and 
%                   Time (1xn). This is the validation data
%
%OUTPUT
%    -  SA_diff   : Matrix of absolute salinity difference between the
%                   model and the instruments. One row for each instrument.
%                   Time in column. Positive means model overshoot   [g/kg]
%    -  SA_mooring: Matrix of absolute salinity data from mooring    [g/kg]
%    -  Y         : Depth of the instruments.                           [m]
%    -  stats     : Table with stats for the whole time.

%VERSION 1.1, Sept 2019
%% DISPLAY
%%disp('Absolute Salinity evaluation')

%% SRUCTURE INFO
SA_e=SA_eval.SAmatrix;
Y_e=SA_eval.Y;
Time_e=SA_eval.Time;

%% Verifications
% Check number of input
if ~(nargin == 8)   %Checking number of args
    error('SA_Eval_fast :  Requires 8 inputs')
end

% Check time inputs
if ~(isnumeric(t0)&&isnumeric(tf))
    error('SA_Eval_fast : t0 and tf have to be in serial time')
end
if mod(((tf-t0)*24*60),dt)~=0
    error('SA_Eval_fast : (tf-t0)/dt has a remainder')
end
if Time_e(1)>t0
    error('SA_Eval_fast : Time_m first stamp too late')
end
if Time_e(end)<tf
    error('SA_Eval_fast : Time_m last stamp too early')
end

% Check spacial input
if mod((y_bot-y_top),dy)~=0
    error('SA_Eval_fast : (y_bot-y_top)/dt has a remainder')
end

% Check Tmatrix dimensions
[m,n]=size(SAmatrix);
if m~=((y_bot-y_top)/dy+1)
    error('SA_Eval_fast : Tmatrix not the right size')
end
if n~=((tf-t0)*24*60/dt)+1
    error('SA_Eval_fast : Tmatrix not the right size')
end 

%% DATA alignement
%time alignement with model
SA_e=interp1(Time_e,SA_e',(t0:dt/24/60:tf));
Y_e=interp1(Time_e,Y_e',(t0:dt/24/60:tf));
SA_e=SA_e';
Y_e=Y_e';
SAmatrix_e=nan(size(SA_e));
%
%vertical alignement to instruments
y=(y_top:dy:y_bot)';   %model depth matrix
for i=1:length(t0:dt/24/60:tf)
    SAmatrix_e(:,i)=interp1(y,SAmatrix(:,i),Y_e(:,i));
end

%outputs
SA_mooring=SA_e;
Y=Y_e;

%% EVALUATION
% Difference
SA_diff=SAmatrix_e-SA_e;
[m,n]=size(SA_diff);
% Stats
if sum(isnan(SA_diff),'all')
    disp('SA_Eval_fast : !WARNING!, NaN in SAdiff')
    if sum(isnan(SA_diff),'all')==length(reshape(SA_diff,1,[]))
        SA_diff(:)=0;  %zero difference if no salinity data...
    end
end
    
Mean=mean(SA_diff,'all','omitnan');
Median=median(SA_diff,'all','omitnan');
STD=std(SA_diff,0,'all','omitnan');
Max=max(SA_diff,[],'all','omitnan');
Min=min(SA_diff,[],'all','omitnan');
RMSE=rms(reshape(SA_diff,1,m*n),'omitnan');
MAE=mean(abs(SA_diff),'all','omitnan');
SStot=sum((SA_e-mean(SA_e,'all','omitnan')).^2,'all','omitnan');
SSres=sum((SA_diff).^2,'all','omitnan');
R2=1-SSres/SStot;
Values=[Mean; Median; STD; Max; Min; RMSE; MAE; R2];

stats=table(Values);
stats.Properties.RowNames={'Mean';'Median';'STD';'Max';'Min';'RMSE';'MAE';'R2'};

end