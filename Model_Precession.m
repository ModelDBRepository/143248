%
% Code to produce main figures in:
%
% "Phase precession through acceleration of local theta rhythm: 
% a biophysical model for the interaction between place cells and local 
% inhibitory neurons", JOURNAL OF COMPUTATIONAL NEUROSCIENCE.
% doi: http://dx.doi.org/10.1007/s10827-011-0378-0
%
% Core model parameters are set in Setup_Parameters
%
% Authors:
% Luisa Castro and Paulo Aguiar
% Faculdade de Ciencias da Universidade do Porto
% contact email: pauloaguiar@fc.up.pt
%
%%

% clean up
clear all
close all

% data containers
FIRINGS_OUT_i=[];
FIRINGS_OUT_e=[];
ISIs_i=[];
ISIs_e=[];
rate_in=[];
rate_out=[];
PHASE_TOT_e=[];
PHASE_TOT_i=[];
vrun_e=[];
vrun_i=[];
vplot_e=[];

nrun=50;     % number of runs

tic
for runy=1:nrun
    
    Setup_Parameters % parameters are set here
    
    FIN(runy)=fini;
    
       for i=2:bins
        t=(i-1)*dt;
        
        
        % Inhibitory neuron receiving the input train
        % Computing current
        g_ti_rise(i)=0;
        g_ti_decay(i)=g_ti_decay(i-1)*exp(-dt/tausin_decay_ti)+FIRINGS_EC(i)*factor*gmax_ti;
        g_ti(i)=g_ti_decay(i)-g_ti_rise(i);
        Isyn_ti(i)=g_ti(i)*(V_i(i-1)-Er_ti);
        
        
        % Computing inhibitory membrane potential
        V_i(i)=V_i(i-1)+(dt/Taum_i)*(Vreset_i-V_i(i-1)+Rm_i*PM(i)-Rm_i*Isyn_ti(i));
        
        
        if V_i(i)>=PARF(1,4)
            firingtime_i=[firingtime_i t];
            vrun_i=[vrun_i runy];
            V_i(i)=PARF(1,3);
            if i+delay/dt<bins
                FIRINGS_INHIB(i+delay/dt)=1;
                
            end
        end
        
        % Excitatory neuron receiving the input train
        % Computing current
        g_te(i)=g_te(i-1)*exp(-dt/tausin_te)+FIRINGS_EC(i)*gmax_te;
        Isyn_te(i)=g_te(i)*(V_e(i-1)-Er_te);
        
        % Excitatory neuron receiving the inhibitory spikes
        % Computing current
        g_ie(i)=g_ie(i-1)*exp(-dt/tausin_ie)+FIRINGS_INHIB(i)*gmax_ie;
        Isyn_ie(i)=g_ie(i)*(V_e(i-1)-Er_ie);
        
        
        % Computing excitatory membrane potential
        V_e(i)=V_e(i-1)+(dt/Taum_e)*(Vreset_e-V_e(i-1)-Rm_e*Isyn_te(i)-Rm_e*Isyn_ie(i));
        
        if V_e(i)>=PARF(2,4)
            firingtime_e=[firingtime_e t];
            vrun_e=[vrun_e runy];
            V_e(i)=PARF(2,3);
        end
        
        ISIunit_i=firingtime_i(2:length(firingtime_i))-firingtime_i(1:length(firingtime_i)-1);
        ISIunit_e=firingtime_e(2:length(firingtime_e))-firingtime_e(1:length(firingtime_e)-1);
        
        
    end
    
    % Computing Phases of the excitatory neuron
    PerTheta=100;                   %[ms]
    phase=mod((2*pi/PerTheta)*mod(firingtime_e,PerTheta)'+fini,2*pi);
    
    % Seting the first spikes to phase 2pi in each run for the excitatory neuron
    if ~isempty(phase)
        vector=2*pi-phase(1);
        phasef=mod(phase+vector,2*pi);
        phasef(find(phasef==0))=2*pi;
        PHASE_TOT_e=[PHASE_TOT_e ;phasef];
    end
    
    
    % Computing Phases of the inhibitory neuron
    phasei=mod((2*pi/PerTheta)*mod(firingtime_i,PerTheta)'+fini,2*pi);
    
    % Setting the first spikes to phase 2pi in each run for the inhibitory neuron
    if ~isempty(phasei)
        vectori=2*pi-phasei(1);
        phasefi=mod(phasei+vectori,2*pi);
        phasefi(find(phasefi==0))=2*pi;
        PHASE_TOT_i=[PHASE_TOT_i ;phasefi];
    end
    
    FIRINGS_OUT_i=[FIRINGS_OUT_i firingtime_i];
    FIRINGS_OUT_e=[FIRINGS_OUT_e firingtime_e];
    ISIs_i=[ISIs_i ISIunit_i];
    ISIs_e=[ISIs_e ISIunit_e];
    
    clear firingtime_e
    clear firingtime_i
    clear phase
    clear phasei
    
end
toc

% Some interesting results
Firing_rate_Inhib=length(FIRINGS_OUT_i)/nrun/T;
Firing_rate_Exc=length(FIRINGS_OUT_e)/nrun/T;
mean_ISIs_inhib=mean(ISIs_i);
std_ISIs_inhib=std(ISIs_i);
mean_ISIs_excit=mean(ISIs_e);
std_ISIs_excit=std(ISIs_e);

%% Plotting fig 3a (or 9 by changing profile in the function intensity_func)
figure
fini_u=0;
x=[0:T];
plot(x,intensity_func(x,profile,fmin,fmax,B_mod, A_mod,fTheta/1000,miu,sig,fini_u,T),'k','Linewidth',2)
xlabel('Time [ms]')
ylabel('Intensity Function')


%% Plotting figs 3b and 3c
figure
subplot(211)
plot(times_sim, PM,'k--',times_sim, -Isyn_ti,'g-', times_sim,-Isyn_te, 'k',times_sim,-Isyn_ie, 'k-.')
xlabel('Time [ms]','fontsize',8)
ylabel('Current [nA]','fontsize',8)
h = legend('PM','I_{synci}','I_{synce}','I_{synie} ',4);
set(h,'Box', 'off','fontsize',8,'Orientation','horizontal','Location','South')
subplot(212)
plot(times_sim, V_i,'g',times_sim, V_e,'k');
xlabel('Time [ms]','fontsize',8)
ylabel('Membrane Potential [mV]','fontsize',8)
ylim([-80 -40])
h = legend('V_{inhibitory}','V_{excitatory}',2);
set(h,'Box','off','fontsize',8,'Orientation','horizontal','Location','North')
hold off

%% Plotting fig 3d
figure
FIRINGS_last=FIRINGS_OUT_e(find(vrun_e==nrun));
Theta=sin(2*pi*fTheta/1000.*times_sim+pi/2);
ThetaFiring=sin(2*pi*fTheta/1000.*FIRINGS_last+pi/2);
plot(times_sim,Theta,'k','Linewidth',2);xlim([0 T+dt]);ylim([-3 3]);hold on
plot(FIRINGS_last,ThetaFiring,'go');
ylabel('Phase Precession - Principal neuron')
xlabel('Time[ms]')
ylim([-1.5 1.5])


%% Plotting fig 4a (or 9 by changing profile in the function intensity_func)
figure
plot(FIRINGS_OUT_e,PHASE_TOT_e,'k.');
ylim([0 2*pi+.3])
xlim([0 T+dt])
set(gca,'YTick',0:pi/2:2*pi);
set(gca,'YTickLabel',{'0','pi/2','pi','3pi/2','2pi'})
xlabel('Time [ms]')
ylabel('Phase [rad] - Exc')


%% Plotting fig 5b
figure
plot(FIRINGS_OUT_i,PHASE_TOT_i,'k.');
ylim([0 2*pi+.3])
xlim([0 T+dt])
set(gca,'YTick',0:pi/2:2*pi);
set(gca,'YTickLabel',{'0','pi/2','pi','3pi/2','2pi'})
xlabel('Time [ms]')
ylabel('Phase [rad] - Inh')

%% Plotting fig 5d
figure
FIRINGS_last=FIRINGS_OUT_i(find(vrun_i==nrun));
Theta=sin(2*pi*fTheta/1000.*times_sim+pi/2);
ThetaFiring=sin(2*pi*fTheta/1000.*FIRINGS_last+pi/2);
plot(times_sim,Theta,'k','Linewidth',2);xlim([0 T+dt]);ylim([-3 3]);hold on
plot(FIRINGS_last,ThetaFiring,'ko');
ylabel('Phase Precession - Local interneuron')
xlabel('Time [ms]')
ylim([-1.2 1.2])


%% Plotting fig 6
figure
hist(FIRINGS_OUT_e/PerTheta,200)
xlabel('Theta Cycles')

% Call additional script to produce other figures
neuron = 'exc';
Firing_Freq;
neuron = 'inh';
Firing_Freq;
