% Set the parameters for the main mfile: Model_Precession.m

T=1000;                             % simulation time duration [ms](T must be >=1000ms and multiple of 100ms)
dt=0.01;                            % time step [ms]
n_int=T/dt;                         
if ~isinteger(n_int)                % ensure that T is multiple of dt
    dt=T/round(n_int);
end
    
times_sim=[0: dt:T-dt];             % simulation times [ms]
bins=length(times_sim);             % number of steps in the simulation


N=1000;                             % number of input units
fmin_avg=0.5;                       % minimum average rate on each input [Hz]
fmax_avg=4.5;                       % maximum average rate on each input [Hz]
fmin=(fmin_avg/1000)*N;             % baseline level for the inputs [kHz]
fmax=(fmax_avg/1000)*N;             % topline level for the inputs [kHz]
B_mod=(fmax-fmin)/fmin;             % scaling factor for the gaussian wave

% Chose the shape of the input stimulus
profile='G';                        % write 'G' for Gaussian; 'S' for Square; 'P' for Pyramid and 'R' for Ramping
sig=T/4;                            % standard deviation value for the gaussian wave [ms]
miu=T/2;                            % mean value for the gaussian wave [ms]

% Inicialize variables
FIN=zeros(1,nrun);
g_ti_rise=zeros(1,bins);
g_ti_decay=zeros(1,bins);
g_ti=zeros(1,bins);
Isyn_ti=zeros(1,bins);
V_i=zeros(1,bins);
g_te=zeros(1,bins);
Isyn_te=zeros(1,bins);
g_ie=zeros(1,bins);
Isyn_ie=zeros(1,bins);
V_e=zeros(1,bins);
FIRINGS_INHIB=zeros(1,bins);

%Aleatory number to set the theta phase in the begining of each run
fini=rand*2*pi;                     % [rad]

%Neurons Parameters
%Fix    (RM,         Taum,   Vreset,   Vth,  Vinicial,   ER)
%       megaohm,     ms      mV        mV     mV        mV
PARF=[200  200  -70 -50  -70  -80 ;200  200  -70 -50  -70  0];
firingtime_i=[];      firingtime_e=[];  % firing times of neurons [ms]
V_i(1)=PARF(1,5);    V_e(1)=PARF(2,5);
Vreset_i=PARF(1,3);  Vreset_e=PARF(2,3);
Rm_i=PARF(1,1);      Rm_e=PARF(2,1);
Taum_i=PARF(1,2);    Taum_e=PARF(2,2);

%Synapses Parameters: te, ti e ie (ie, train->inhib; train->excit; inhib->excit)
tausin_rise_ti=2;       % [ms]
tausin_decay_ti=5;      % [ms]
g_ti_rise(1)=0;         % [microS]
g_ti_decay(1)=0;        % [microS]
g_ti(1)=0;              % [microS]
Isyn_ti(1)=0;           % [nA]
gmax_ti=0.03/1000;      % [microS]
Er_ti=0;                % [mV]
factor=1;

tausin_te=5;            % [ms]   
g_te(1)=0;              % [microS]
Isyn_te(1)=0;           % [nA]
gmax_te=.4/1000;        % [microS]
Er_te=0;                % [mV]

tausin_ie=20;           % [ms]
g_ie(1)=0;              % [microS]
Isyn_ie(1)=0;           % [nA]
gmax_ie=.04;            % [microS]
Er_ie=PARF(1,6);        % [mV]
delay=5;                % [ms]

%Pacemaker
fTheta=10;              % Theta frequency [Hz]
k=0.1;                  % [nA]
base=0.25;              % [nA]
PM=k*cos(2*pi*((fTheta)/1000)*times_sim+fini)+base;
A_mod=0.2;             % sinusoidal suavization factor - value in [0,1]
% use A_mod=0 to simulate profile in figure 9a

% Spike train modulated by profile set in intensity_func function
% To produce figs 8a, 8b and 8c set fmin=2.5(eg) and fmax=fmin in the begining so that B_mod=0
[trainin]=PoissonProc_Generator(T,dt,profile,fmin,fmax,B_mod,A_mod,fTheta/1000,miu,sig,fini);

% Trainin produced above is a vector with fraction numbers in [0,T]
% (4 decimal places)
% The next lines turn it into a vector (FIRINGS_EC) of size T/dt, with integer
% entrances meaning the number of spikes occurring in that milisecond.
% This will be the input in the integrate and fire equation in Model_Precession.m

AUX1=round(trainin/dt);
AUX2=zeros(1,T/dt+1);
for i=1:length(AUX1)
    AUX2(AUX1(i)+1)=AUX2(AUX1(i)+1)+1;
end
FIRINGS_EC=AUX2;

