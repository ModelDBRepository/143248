% Run this m.file after Model_Precession.m to produce other figures, for the
% excitatory neuron or for the inhibitory neuron
% Luisa Castro, FCUP

% which figures should be produced?
if exist('neuron') == 0
    % should be either 'exc', for the excitatory neuron, or 'inh' for the inhibitory neuron
    neuron = 'inh';
end

clear firestr
clear ge
clear gig
clear mean_t
clear mean_p
clear std_t
clear std_p

if neuron == 'exc'
    aa=FIRINGS_OUT_e;
    bb=PHASE_TOT_e;
    cc=vrun_e;
    ncl=7;              % number of clusters visualized in the phase precession plot (fig 4a)
    dd=0;
elseif neuron == 'inh'
    PHASE_TOT_ir=unwindphases(FIRINGS_OUT_i,PHASE_TOT_i,T/2,pi);    % set for plot figs 5a and 5c
    aa=FIRINGS_OUT_i;
    bb=PHASE_TOT_ir;
    cc=vrun_i;
    ncl=11;              % number of clusters visualized in the phase precession plot (fig 5b)
    dd=-pi/4;
end


% Saving data in structures
firestr = struct('firet',[] ,'ISI',[],'ptomed',[],'firert',[],'phase',[]);
for i=1:max(cc)
    firestr(i).firet=aa(find(cc==i));
    firestr(i).phase=bb(find(cc==i));
end

% Extracting Inter Spike Intervals (ISIs)
for i=1:max(cc)
    gu=length(firestr(i).firet);
    firestr(i).ISI=firestr(i).firet(2:gu)-firestr(i).firet(1:gu-1);
end

% Computing mean points for each consecutive pair of spikes
for i=1:max(cc)
    gi=length(firestr(i).firet);
    firestr(i).ptomed=firestr(i).firet(1:gi-1)+firestr(i).ISI/2;
end


% Computing firing rates by inverting ISI
for i=1:max(cc)
    firestr(i).firert=1000./firestr(i).ISI;
end

% Computing the max number of spikes of the simulations
for i=1:max(cc)
    ge(i)=length(firestr(i).firet);
end
ge=max(ge);

% Computing the mean and standard deviation for ptomed, firert, phase and
% firet from all the runs
valxm=[];
valxstd=[];
valym=[];
valystd=[];
valzm=[];
valzstd=[];
valtm=[];
valtstd=[];


for j=1:ge
    auxx=[];
    auxy=[];
    auxz=[];
    auxt=[];
    for i=1:max(cc)
        if length(firestr(i).ptomed)>=j
            auxx=[auxx firestr(i).ptomed(j)];
            auxy=[auxy firestr(i).firert(j)];
        end
        if length(firestr(i).firet)>=j
            auxz=[auxz firestr(i).phase(j)];
            auxt=[auxt firestr(i).firet(j)];
        end
    end
    valxm(j)=mean(auxx);
    valxstd(j)=std(auxx);
    valym(j)=mean(auxy);
    valystd(j)=std(auxy);
    valzm(j)=mean(auxz);
    valzstd(j)=std(auxz);
    valtm(j)=mean(auxt);
    valtstd(j)=std(auxt);
end

valxm=valxm(find(~isnan(valxm)));
valxstd=valxstd(find(~isnan(valxstd)));
valym=valym(find(~isnan(valym)));
valystd=valystd(find(~isnan(valystd)));
valzm=valzm(find(~isnan(valzm)));
valzstd=valzstd(find(~isnan(valzstd)));
valtm=valtm(find(~isnan(valtm)));
valtstd=valtstd(find(~isnan(valtstd)));
tmax=max(aa);
tmin=min(aa);

% Comment from here to run Fic_Aux_TxInhib.m

if neuron == 'exc'
    % Ploting fig 4c
    figure
    errorbar([tmin valxm tmax],[0 valym 0],[0 valystd 0],'black','LineWidth',2)
    xlim([0 T+dt])
    xlabel('Time [ms]')
    ylim([0 15])
    ylabel(['Firing Rate [Hz]', ' - ', neuron])
    hold on
    herrorbar([tmin valxm tmax],[0 valym 0],[0 valxstd 0],'k')
    hold off
    
elseif neuron == 'inh'
    % Ploting fig 5a
    figure
    errorbar([0 valxm T],[10 valym 10],[0.25 valystd 0.25],'black','LineWidth',2)
    xlim([0 T+dt])
    xlabel('Time [ms]')
    ylim([0 15])
    ylabel(['Firing Rate [Hz]', ' - ', neuron])
    hold on
    herrorbar([0 valxm T],[10 valym 10],[0 valxstd 0],'k')
    hold off
end


% Defining clusters to plot fig 4b
dispas=aa;
phases=bb;
IDX = clusterdata(dispas','maxclust',ncl);
clust=struct('cl_firet',[] ,'cl_phase',[]);
for i=1:ncl
    clust(i).cl_firet=dispas(find(IDX==i));
    clust(i).cl_phase=phases(find(IDX==i));
end

% Plotting clusters of firing phases
figure
ylim([0 2*pi+.3])
xlim([0 T+dt])
plot(clust(1).cl_firet,clust(1).cl_phase,'g.'); hold on
plot(clust(2).cl_firet,clust(2).cl_phase,'b.'); hold on
plot(clust(3).cl_firet,clust(3).cl_phase,'c.'); hold on
plot(clust(4).cl_firet,clust(4).cl_phase,'m.'); hold on
plot(clust(5).cl_firet,clust(5).cl_phase,'y.'); hold on
plot(clust(6).cl_firet,clust(6).cl_phase,'r.'); hold on
plot(clust(7).cl_firet,clust(7).cl_phase,'k.'); hold on
if ncl>7
    hold on
    plot(clust(8).cl_firet,clust(8).cl_phase,'g.');
    plot(clust(9).cl_firet,clust(9).cl_phase,'b.'); hold on
    plot(clust(10).cl_firet,clust(10).cl_phase,'c.'); hold on
    plot(clust(11).cl_firet,clust(11).cl_phase,'m.');
end
ylim([0 2*pi+.3])
xlim([0 T+dt])
set(gca,'YTick',0:pi/2:2*pi);
set(gca,'YTickLabel',{'0','pi/2','pi','3pi/2','2pi'})
xlabel('Time [ms]')
ylabel(['Phase [rad]', ' - ', neuron])
hold off

% Computing the mean and standard deviation values for firing times and
% phases in each cluster
for i=1:ncl
    mean_t(i)=mean(clust(i).cl_firet);
    mean_p(i)=mean(clust(i).cl_phase);
    std_t(i)=std(clust(i).cl_firet);
    std_p(i)=std(clust(i).cl_phase);
end

gig=[mean_t; std_t;mean_p; std_p]';
gig=sortrows(gig);                      % sort by spike times

% Plotting fig 4b (exc) or 5c (inh)
figure
errorbar(gig(:,1) , gig(:,3), gig(:,4) ,'black','LineWidth',2)
ylim([dd 2*pi+.3])
xlim([0 T+dt])
set(gca,'YTick',0:pi/2:2*pi);
set(gca,'YTickLabel',{'0','pi/2','pi','3pi/2','2pi'})
xlabel('Time [ms]')
ylabel(['Firing Phase [rad]', ' - ', neuron])
hold on
herrorbar(gig(:,1) , gig(:,3), gig(:,2) ,'k')
hold off


