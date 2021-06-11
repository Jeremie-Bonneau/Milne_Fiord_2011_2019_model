%% Milne Fiord Epishelf Lake model for 2018-2019
% see README file for more details
% do not forget to add the folder's path

% change the directory to point at the right place
cd 'C:\Users\jerem\OneDrive\Documents\UBC\archive\Milne_Fiord_Epishelf_Lake_Model'
%% VALUES
t0=datenum(2018,09,05);  % start time
tf=datenum(2019,05,29);  % final time
dt=24*60;                % 1 day time steps
y_top=3.5;               % top of grid, m
y_bot=25;                % bottom of grid, m
dy=0.1;                  % grid space, m
time=t0:dt/24/60:tf;     % time array
y=[y_top:dy:y_bot]';

%% EVALUATION DATA SETS
% temperature for CT evaluation
Tmatrix_eval=importdata('data\Tmatrix_0_30_day_10cm.mat');
% some trim to make it faster
Teval.Tmatrix=interp1(Tmatrix_eval.Time,Tmatrix_eval.Tmatrix',t0-1:dt/60/24:tf+1);
Teval.Tmatrix=Teval.Tmatrix';
Teval.Time=t0-1:dt/24/60:tf+1;
Teval.Y=Tmatrix_eval.Y;
% SA daily data for evaluation
SA_data=importdata('data\SA_day.mat');
% some trim to make it faster
SAeval.SAmatrix=SA_data.temp_18_19;
SAeval.Y=SA_data.Y_18_19;
SAeval.Time=SA_data.time_18_19;

%% PREPROCESS
% Mesh
empty_mesh=meshing(t0,tf,dt,y_top,y_bot,dy);
% Initial consitions
[SAmatrix,CTmatrix] = IC_mooring_ctd(empty_mesh,t0,y_top,y_bot,dy);
% Boundary conditions
SAmatrix(1,2:end)=0;                         %salt flux at top = 0
SAmatrix(end,2:end)=0;                       %salt flux at bottom = 0
CTmatrix=TBC(CTmatrix,t0,tf,dt,y_top,y_bot); %kowned temperature at top and bottom

%% PROCESS
% The model is run with with a pair of outflow coefficient and the best
% daily mixing coefficients are found for everyday. Then, the whole model
% timeseries is evaluated and the coefficient of agreement is written in
% Ca. Then another year is simulated using a different pair of mixing coef.
% The model output with the best coefficient of agreement is written in a
% structure. 

% Parameters values
h0=[7.8];                   %minimum ice shelf draft array
coef=[5.5];                                  %friction*width coef array

h0matrix=repmat(h0,length(coef),1);
coefmatrix=repmat(coef',1,length(h0));

% Some variables
SA=cell(1,length(h0)*length(coef));        % cell to store all SAmatrix from find_daily_K
CT=cell(1,length(h0)*length(coef));        % cell to store all CTmatrix from find_daily_K
K_CT_top=cell(1,length(h0)*length(coef));  % cell to store all top heat eddy coeffcients from find_daily_K
K_CT_bot=cell(1,length(h0)*length(coef));  % cell to store all bottom heat eddy coeffcients from find_daily_K
K_SA_top=cell(1,length(h0)*length(coef));  % cell to store all top salt eddy coeffcients from find_daily_K
K_SA_bot=cell(1,length(h0)*length(coef));  % cell to store all bottom salt eddy coeffcients from find_daily_K
Outflow=cell(1,length(h0)*length(coef));   % cell to store all outflow values from find_daily_K
Ca=nan(1,length(h0)*length(coef));         % array to store Ca

tic
% Computation for all coef, h0, K_t_top and K_t_bot
parfor i=1:(length(h0)*length(coef))
    [SA{i},CT{i},K_CT_top{i},K_CT_bot{i},K_SA_top{i},K_SA_bot{i},Outflow{i}] = find_daily_K(t0,tf,y_top,y_bot,dy,SAmatrix,CTmatrix,coefmatrix(i),h0matrix(i),SAeval,Teval);
end
% Evaluation of the best parameters
for i=1:length(SA)
    Ca(i)=model_score(t0,tf,dt,y_top,y_bot,dy,SA{i},CT{i},SAeval,Teval);
end

%best fitting parameters
[~,idx]=max(Ca);
SAmatrix=SA{idx};
CTmatrix=CT{idx};
K_CT_top=K_CT_top{idx};
K_CT_bot=K_CT_bot{idx};
K_SA_top=K_SA_top{idx};
K_SA_bot=K_SA_bot{idx};
Outflow=Outflow{idx};

Ca=reshape(Ca,length(coef),length(h0));

toc
%% VISUAL
figure; 
contourf(CTmatrix,20);
colormap jet
axis ij

figure;
contourf(SAmatrix,20);
colormap jet
axis ij

%% 
% figure;
% contourf(h0,coef,Ca);

%% POSTPROCESS
[CTdiff,CTmooring,CTstats] = CT_Eval_fast(t0,tf,24*60,y_top,y_bot,dy,SAmatrix,CTmatrix,Teval);
[SAdiff,SAmooring,Ymooring,SAstats] = SA_Eval_fast(t0,tf,24*60,y_top,y_bot,dy,SAmatrix,SAeval);

%% structure
weir_05sept2018_29may2019.SAmatrix=SAmatrix;
weir_05sept2018_29may2019.SAdiff=SAdiff;
weir_05sept2018_29may2019.SAmooring=SAmooring;
weir_05sept2018_29may2019.SAstats=SAstats;
weir_05sept2018_29may2019.SAmooring_Y=Ymooring;
weir_05sept2018_29may2019.CTmatrix=CTmatrix;
weir_05sept2018_29may2019.CTdiff=CTdiff;
weir_05sept2018_29may2019.CTmooring=CTmooring;
weir_05sept2018_29may2019.CTstats=CTstats;
weir_05sept2018_29may2019.outflow=Outflow;
weir_05sept2018_29may2019.time=time;
weir_05sept2018_29may2019.K_CT_top=K_CT_top;
weir_05sept2018_29may2019.K_CT_bot=K_CT_bot;
weir_05sept2018_29may2019.K_SA_top=K_SA_top;
weir_05sept2018_29may2019.K_SA_bot=K_SA_bot;
weir_05sept2018_29may2019.Ca=Ca;
weir_05sept2018_29may2019.h0=h0;
weir_05sept2018_29may2019.coef=coef;
save 05sept2018_01jun2019_h078_Cb55.mat weir_05sept2018_29may2019