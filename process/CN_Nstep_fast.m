function [Profile_out] = CN_Nstep_fast(Profile_in,K,BC_top,value_top,BC_bot,value_bot,dt,N,y_top,y_bot,dy)
% CN_Nstep solves N time steps dt of the 1D transport eddy diffusivity 
% equation, no advection term, using the Crank-Nicolson scheme with
% specified top and bottom boundary conditions.
% good source: http://www.nada.kth.se/~jjalap/numme/FDheat.pdf

% Using n for time, i for space (vertical):

% T_i,n - T_i,n-1               T_i-1,n - 2*T_i,n + T_i+1,n
% _______________ = [Kt + Km] * ___________________________ +
%    delta_t                           2*delta_y^2
%
%             T_i-1,n-1 - 2*T_i,n-1 + T_i+1,n-1
% [Kt + Km] * ________________________________
%                      2*delta_y^2
%
%    Kt         : Turbulent diffusivity coefficient
%    Km         : Molecular diffusivity coefficient
%
% INPUTS
%    Profile_in : nx1 array of scalar values at time t-dt
%    K          : nx1 array of the combined diffusivities           [m^2/s]
%    BC_top     : Top boundary condition, 'D'(Diriclet) or 'N'(Neumann)
%    value_top  : Top BC value
%    BC_bot     : Bottom boundary condition, 'D'(Diriclet) or 'N'(Neumann)
%    value_bot  : Bottom BC value
%    dt         : Time discretization                                 [min]
%    N          : Number of time steps
%    y_top      : Top depth                                             [m]
%    y_bot      : Bottom depth                                          [m]
%    dy         : Vertical discretization                               [m]
%
% OUTPUT
%    T_model    : nx1 array of the solved 1D transport equation at time t
%
% VERSION 1, August 2019

%% DISPLAY
%disp('CN_1step : Solving 1 step of 1D transport eq with Crank-Nicolson scheme')

%% VERIFICATIONS
% comment to make it faster
% % Check number of input
% if ~(nargin == 11)   %Checking number of args
%     error('CN_Nstep :  Requires 11 inputs')
% end
% 
% % Check spacial input
% if mod((y_bot-y_top),dy)~=0
%     error('CN_Nstep : (y_bot-y_top)/dy has a remainder')
% end
%  
% % Check input K
% if sum(isnan(K))~=0
%     error('CN_Nstep : error in K: has NaN value(s)')
% end
% 
% % Check K_matrix
% if sum(isnan(Profile_in),'all')~=0
%     error('CN_Nstep : error in Profile_in: has NaN value(s)')
% end
% if size(K)~=size(Profile_in)
%     error('CN_Nstep : Profile_in and K not same size')
% end
% 
% % Check BC
% if ~(BC_top=='D'||BC_top=='N')
%     error('CN_Nstep : BC_top format not recognized')
% end
% if ~(BC_bot=='D'||BC_bot=='N')
%     error('CN_Nstep : BC_top format not recognized')
% end

%% Matrix system
% Solving matrix system for one time step:
%     Coef_matrix*T_n = a*T_i-1,n-1 + (1/dt+2a)*T_i,n-1 + a*T_i+1,n1
%                 A*x = B

nx=length(Profile_in);        % number of vertical nodes
dt=dt*60;                     % transfering in seconds


%% 
X=Profile_in;
for i=1:N
    %% Building the coefficients matrix A
    % main diagonal
    diag=(K./(dy^2))+1/dt;             % diagonal coefficients column
    diag=eye(nx).*diag;                % matrix with diag coefficients
    % side-diagonal
    s_diag=(K./(-dy^2)./2);            % side-diag coefficients column
    index=zeros(nx,1); index(2)=1;     % row to build matrix with 1 at the side diag positions
    s_diag=toeplitz(index).*s_diag;    % matrix with the 2 side diagonals
    % total
    A=diag+s_diag;
    %
    %% Building vector B
    % B = a*T_i-1,n-1 + (1/dt-2a)*T_i,n-1 + a*T_i+1,n
    a=K./(dy^2)./2;                % coefficient a
    B1=[0; X(1:end-1)].*a;         % a*T_i-1,n-1
    B2=X.*(-2*a+1/dt);             % (1/dt+2a)*T_i,n-1
    B3=[X(2:end); 0].*a;           % a*T_i+1,n-1
    B=B1+B2+B3; 

    %% Boundary conditions
    %
    % Diriclet
    % T_i,n = Value given
    %
    % Neumann
    % second order accuracy derivative
    %
    %         -3T_i,n + 4T_i+1,n - 1T_i+2,n
    % dT_dy = _____________________________ = value given
    %                     2dy
    % 
    if strcmp('D',BC_top)
        A(1,2)=0;  A(1,1)=1;
        B(1)=value_top;
    end
    if strcmp('D',BC_bot)
        A(end,end-1)=0; A(end,end)=1;
        B(end)=value_bot;
    end
    if strcmp('N',BC_top)
        A(1,1)=-3; A(1,2)=4; A(1,3)=-1;
        B(1)=value_top*2*dy;
    end
    if strcmp('N',BC_bot)
        A(end,end)=-3; A(end,end-1)=4; A(end,end-2)=-1; 
        B(end)=value_bot*2*dy;
    end

    %% Solving system
    X=A\B;
end

Profile_out=X;
    
end