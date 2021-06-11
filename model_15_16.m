%% Milne Fiord Epishelf Lake model for 2015-2016
% see README file for more details
% do not forget to add the folder's path

% change the directory to point at the right place
cd 'C:\Users\jerem\OneDrive\Documents\UBC\archive\Milne_Fiord_Epishelf_Lake_Model'
%% VALUES
t0=datenum(2015,09,06); % start time
tf=datenum(2016,06,04); % final time
dt=24*60;               % 1 day time steps
y_top=3.0;                % top of grid, m
y_bot=25;               % bottom of grid, m
dy=0.1;                 % grid space, m
time=t0:dt/24/60:tf;    % time array
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
SAeval.SAmatrix=SA_data.temp_15_16;
SAeval.Y=SA_data.Y_15_16;
SAeval.Time=SA_data.time_15_16;

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
% For each day, the model is run using all the combination of parameters
% values and an evaluation score is given. Once all the parameters
% combination have been tried out, the one with the best score (best fit)
% is use to solve for the day and written in CTmatrix and SAmatrix. Then on
% to the next day. 

% Parameters values
h0=[5.3];                   %minimum ice shelf draft
coef=[3.4];                                  %friction*width coef

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

% figure;
% contourf(h0,coef,Ca);

%% POSTPROCESS
[CTdiff,CTmooring,CTstats] = CT_Eval_fast(t0,tf,24*60,y_top,y_bot,dy,SAmatrix,CTmatrix,Teval);
[SAdiff,SAmooring,Ymooring,SAstats] = SA_Eval_fast(t0,tf,24*60,y_top,y_bot,dy,SAmatrix,SAeval);

%% structure
weir_06sept2015_04jun2016.SAmatrix=SAmatrix;
weir_06sept2015_04jun2016.SAdiff=SAdiff;
weir_06sept2015_04jun2016.SAmooring=SAmooring;
weir_06sept2015_04jun2016.SAstats=SAstats;
weir_06sept2015_04jun2016.SAmooring_Y=Ymooring;
weir_06sept2015_04jun2016.CTmatrix=CTmatrix;
weir_06sept2015_04jun2016.CTdiff=CTdiff;
weir_06sept2015_04jun2016.CTmooring=CTmooring;
weir_06sept2015_04jun2016.CTstats=CTstats;
weir_06sept2015_04jun2016.outflow=Outflow;
weir_06sept2015_04jun2016.time=time;
weir_06sept2015_04jun2016.K_CT_top=K_CT_top;
weir_06sept2015_04jun2016.K_CT_bot=K_CT_bot;
weir_06sept2015_04jun2016.K_SA_top=K_SA_top;
weir_06sept2015_04jun2016.K_SA_bot=K_SA_bot;
weir_06sept2015_04jun2016.Ca=Ca;
weir_06sept2015_04jun2016.h0=h0;
weir_06sept2015_04jun2016.coef=coef;
save 06sept2015_04jun2016_h053_Cb34.mat weir_06sept2015_04jun2016