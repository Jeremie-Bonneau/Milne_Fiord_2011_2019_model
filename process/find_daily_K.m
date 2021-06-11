function [SAmatrix,CTmatrix,K_CT_top,K_CT_bot,K_SA_top,K_SA_bot,Out] = find_daily_K(t0,tf,y_top,y_bot,dy,SAmatrix,CTmatrix,coef,h0,SAeval,Teval)
% function returning SAmatrix and CTmatrix for the specified time interval
% looping through the possible values of mixing coefficients (inside
% function). This function uses the inverse weir flow equation
%
%INPUTS
%    t0            : Initial time                                  [serial]
%    tf            : Final time                                    [serial]
%    y_top         : top of the water column under sea surface          [m]
%    y_bot         : bottom of the water column under sea surface       [m]
%    dy            : Vertical discretization                            [m] 
%    SAmatrix      : SA matrix with IC and BC                        [g/kg]
%    CTmatrix      : CT matrix with IC and BC                        [degC]
%    coef          : friction coefficient * width of channel            [m]
%    h0            : top depth of the outflow                           [m]
%    SAeval        : SA structure for evaluation
%    Teval         : Temperature structure for evaluation
%
%OUTPUT
%    SAmatrix      : SA matrix for the specified time interval       [g/kg]
%    CTmatrix      : CT matrix for the specified time interval       [degC]
%
%VERSION 1, sept 2019
%
%% VERIFICATIONS

%% SOME VALUES
K_CT_top_array=[1.40001 2 4 8 16 32 64 128 256 512 1024].*10^(-7);  %mol+tur heat coef top
K_CT_bot_array=[1.40001 2 4 8 16 32 64 128 256 512 1024].*10^(-7);  %mol+tur heat coef bottom
K_CT_top=nan(1,tf-t0);
K_CT_bot=nan(1,tf-t0);
K_SA_top=nan(1,tf-t0);
K_SA_bot=nan(1,tf-t0);
Out=nan(1,tf-t0);

%% LOOPS

for t=1:tf-t0 %loop through days
    N_2=N2_profile(SAmatrix(:,t),CTmatrix(:,t),y_top,y_bot,dy);
    score=nan(length(K_CT_top_array),length(K_CT_bot_array));     %2d matrix of evaluation score
    [SA,CT,Out(t)] = Outflow_weir(SAmatrix(:,t),CTmatrix(:,t),N_2,y_top,y_bot,dy,coef,h0);
    
            for i=1:length(K_CT_top_array) %loop through possible top turbulent coef   
                for ii=1:length(K_CT_bot_array) %loop through possible bottom turbulent coef
                    % Compute mixing coefficients arrays
                    [Ks,Kh] = K_N2_custom_fast(N_2,y_top,y_bot,dy,K_CT_top_array(i),K_CT_bot_array(ii));
                    % SOLVE 1 DAY
                    [CTday] = CN_Nstep_fast(CT,Kh,'D',CTmatrix(1,t+1),'D',CTmatrix(end,t+1),60,24,y_top,y_bot,dy); %1 day, 30 min time steps
                    [SAday] = CN_Nstep_fast(SA,Ks,'N',SAmatrix(1,t+1),'N',SAmatrix(end,t+1),60,24,y_top,y_bot,dy); %1 day, 30 min time steps
                    % EVALUATION : Coefficient of agreement
                    Ca=model_score(t0+t,t0+t,30,y_top,y_bot,dy,SAday,CTday,SAeval,Teval);
                    score(i,ii)=Ca;
                    if ii>2&&score(i,ii)<score(i,ii-1)
                        break
                    end
                end %ii loop
                if i>2&&score(i,ii)<score(i-1,ii)
                    break
                end
            end %i loop

    % BEST FIT
    % Indexes of the best fit (max coefficient of agreement)
    [~,idx]=max(reshape(score,1,[]));
    [i_idx, ii_idx]=ind2sub(size(score),idx);  %indexes of the best fit
    % Compute mixing coefficients arrays of the best fit
    [Ks,Kh] = K_N2_custom_fast(N_2,y_top,y_bot,dy,K_CT_top_array(i_idx),K_CT_bot_array(ii_idx));
    % Solving with the best fit parameters
    [CTday] = CN_Nstep_fast(CT,Kh,'D',CTmatrix(1,t+1),'D',CTmatrix(end,t+1),30,48,y_top,y_bot,dy); %1 day, 30 min time steps
    [SAday] = CN_Nstep_fast(SA,Ks,'N',SAmatrix(1,t+1),'N',SAmatrix(end,t+1),30,48,y_top,y_bot,dy); %1 day, 30 min time steps
    % Write in output matrices
    SAmatrix(:,t+1)=SAday;
    CTmatrix(:,t+1)=CTday;
    % store parameters value
    K_CT_top(t)=K_CT_top_array(i_idx);
    K_CT_bot(t)=K_CT_bot_array(ii_idx);
    K_SA_top(t)=K_SA_from_K_CT(K_CT_top(t));
    K_SA_bot(t)=K_SA_from_K_CT(K_CT_bot(t));
end %t loop
    
end
