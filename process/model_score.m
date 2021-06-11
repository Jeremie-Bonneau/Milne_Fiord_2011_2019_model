function [Ca] = model_score(t0,tf,dt,y_top,y_bot,dy,SAmatrix,CTmatrix,SAeval,Teval)
% function to compute the coefficient of agreement (Ca) between the model
% and the available data from the mooring. Ca is calculated followingly:
%
%
%                  n_t     RMSE_CT,mod     n_c     RMSE_SA,mod
%         Ca = ( _______ * ___________ + _______ * ___________ )^-1
%                n_t+n_c    STD_CT,moo    n_t+n_c   STD_SA,moo
%
% where
%    -  n_t is the number of temperature measurements on the mooring line
%    -  n_c is the number of conductivity measurements
%    -  RMSE_CT,mod is the RMSE of the model compared to the mooring data
%    -  RMSE_SA,mod is the RMSE of the model compared to the mooring data
%    -  STD_CT,moo is the CT standard deviation of the mooring data
%    -  SDD_SA,moo is the SA standard deviation of the mooring data
%
%INPUT
%    -  t0        : Initial time                                   [serial]
%    -  tf        : Final time                                     [serial]
%    -  dt        : Time discretization                               [min]
%    -  y_top     : Top depth                                           [m]
%    -  y_bot     : Bottom depth                                        [m]
%    -  dy        : Vertical discretization                             [m]
%    -  CTmatrix  : Temperature matrix returned by the model         [degC]
%    -  SAmatrix  : Absolute salinity returned by the model          [g/kg]
%
%OUTPUT
%    -  Ca        : Coefficient of agreement


%% WEIGHTS
if t0>=734641&&tf<=735021
    if t0>=734641&&tf<=734994 %May 2011-May2012
        SAw=5/23;        %n_t/(n_t+n_c)
        CTw=1-SAw;       %n_c/(n_t+n_c)
    elseif t0>734994&&tf<735021 %May 2012 to June 2012
        SAw=3/17;
        CTw=1-SAw; 
    else
        SAw=5/17; %2 or 5 inst (solinsts stopped in Jan)
        CTw=1-SAw;
    end
elseif t0>=735004&&tf<=735054
    SAw=5/23;
    CTw=1-SAw;
elseif t0>=735056&&tf<=735420
    SAw=2/9;
    CTw=1-SAw;
elseif t0>=735421&&tf<=735435
    SAw=3/17;
    CTw=1-SAw;
elseif t0>=735435&&tf<=735793
    SAw=2/11;
    CTw=1-SAw;
elseif t0>=735797&&tf<=736155
    SAw=2/13;
    CTw=1-SAw;
elseif t0>=736158&&tf<=736175
    SAw=4/15;
    CTw=1-SAw;
elseif t0>=736175&&tf<=736527 %2015-2016
    SAw=4/16;
    CTw=1-SAw;
elseif t0>=736533&&tf<=736901
    SAw=5/17;
    CTw=1-SAw;
elseif t0>=736902&&tf<=737249 %2017-2018
    if t0>=737065&&tf<=737249 %jan2018-2018
        SAw=2/17;   %2 inst from jan to end (solinsts stopped in Jan)
        CTw=1-SAw;
    elseif t0>=736902&&tf<737065 %July2017 to Jan2018
        SAw=5/17;   %5 inst from July to jan (solinsts stopped in Jan)
        CTw=1-SAw; 
    else
        SAw=3.5/17; %2 or 5 inst (solinsts stopped in Jan)
        CTw=1-SAw;
    end
elseif t0>=737253&&tf<=737621
    SAw=3/16;
    CTw=1-SAw;
else 
    error('model_score : could not evaluate evaluation weights')
end

%% EVALUATION
[~,~,CTstats] = CT_Eval_fast(t0,tf,dt,y_top,y_bot,dy,SAmatrix,CTmatrix,Teval);
[~,~,~,SAstats] = SA_Eval_fast(t0,tf,dt,y_top,y_bot,dy,SAmatrix,SAeval);

%% Ca
CTscore=CTstats.Values(6)/std(CTmatrix,0,'all').*CTw;     %RMSE normalized by the standard deviation of the model data (i.e. CT data)
SAscore=SAstats.Values(6)/std(SAmatrix,0,'all').*SAw;     %RMSE normalized by the standard deviation of the model data (i.e. SA data)
Ca=1/(CTscore+SAscore);

