%% 0. Prepare paths 
p.scripts = fileparts(fileparts(matlab.desktop.editor.getActiveFilename)); % Store the path of this script
addpath(genpath(p.scripts))

p.pBOSC = fileparts(p.scripts);
p.data = fullfile(p.pBOSC, 'Simresults');
p.files = fullfile(p.pBOSC, 'Files');

% Toolbox paths
addpath('Z:\Toolbox\fieldtrip-20230118') % Path to Fieldtrip
ft_defaults

% Participant and session list
[Nsub, subs, sess] = Omega_subs();

% Frequency parameters
frex = exp(0.6:0.1:3.7); % 1.8 Hz to 40 Hz
Nfrex = length(frex);

source = Omega_source(3423, 1.5);
inside = find(source.inverse.inside);

% Example Simulation Parameters
Nsources =[2924];
Nfreqs =  [10];
Ncycles = [10];
Nsnr = [1];

%% Figure 1. Simulate the signal at source level

cfg = [];
cfg.numgens = 5;
cfg.repetitions = 1;
cfg.Nsources = Nsources;
cfg.seed = 1; % seed for replication
evalc('apgens = Sim_apgenerators(cfg)');

cfg = [];
cfg.dip = [Nsources];
cfg.apgens = apgens;
cfg.simfreq = [Nfreqs];
cfg.simcycles = [Ncycles];
cfg.simsnr = [Nsnr];
cfg.signaldur = 6;
cfg.figflag = 1; % Plot the figures
evalc('simsignal = pBOSC_SimulateSignal(cfg)');

%% Figure 2. Source-reconstruction
evalc('simsignal_source = pBOSC_SimulateBeamformer(simsignal)');

figure('Position',[0 0 1000 800])
plot(source.inverse.pos(inside,2), source.inverse.pos(inside,3),'ok', 'MarkerSize',15, 'LineWidth',2)
set(gca, 'Ylim', [-12 14])
set(gca, 'Xlim', [-14 10])
axis off

figure, plot(simsignal_source.trial{1}(Nsources,:),'k'), axis off,

%% Figure 3. 
cfg = [];
cfg.dip = [Nsources];
cfg.winlen = 2;
cfg.figflag = 0;
evalc('sim_aperiodic = pBOSC_SimulateAperiodic(cfg, simsignal_source)');

%% Figure 4.
cfg = [];
cfg.dip = Nsources;
cfg.ff = Nfreqs;
cfg.cyc = Ncycles;
cfg.snr = Nsnr;
cfg.figflag = 1;
evalc('[model] = pBOSC_SimulatedEpisodes(cfg, simsignal_source, sim_aperiodic)');

%% Figure 5

load('Z:\OMEGA\Enrique\Results\newoscillatory_results.mat')

% Natural Frequencies
all_natfreq = normalize(freq_osc_time_cycles,2,"zscore","robust");
avg_natfreq = squeeze(mean(all_natfreq,1));
avg_domfreq = squeeze(mean(freq_osc_time_cycles,1));
nf=[];
df = [];

% Frequency parameters
frex = exp(0.6:0.1:3.7); % 1.8 Hz to 40 Hz
Nfrex = length(frex);

source = Omega_source(1925,1.5);
voxin = find(source.inside);
Nvoxin = length(voxin);
for v = 1:Nvoxin
    [pks,locs] = findpeaks(avg_natfreq(v,:));
    pks2 = pks(~isinf(pks));
    locs = locs(~isinf(pks));
    temp=frex(locs(find(pks2==max(pks2))));
    nf(v)=temp(1);
end

addpath(genpath('Z:\OMEGA\OMEGA-NaturalFrequencies-main'))
% pBOSC
colcfg.colim   = [0.7 3.4];  
colcfg.colmap   = 'jet_omega_mod';
pBOSC_sourcefig(log(nf)', 'new_natfreq', colcfg)

% Almu
load Z:\OMEGA\Enrique\Results\natfreq_aalMU.mat
natfmu_Almu      = median(natfMU.mu,2,"omitmissing");
colcfg.colim   = [0.7 3.4];  
colcfg.colmap   = 'jet_omega_mod';
pBOSC_sourcefig(log(natfmu_Almu), 'test',colcfg)

load('Z:\OMEGA\Enrique\Results\natfreq_aal_2.mat')
diffs = abs(natfmu_aal - natfmu_aalMU);

diffsnorm = ( diffs-min(diffs) ) / ( max(diffs) - min(diffs));

% Heatmap
natfmu_aal = round(natfmu_aal,1);
natfmu_aalMU = round(natfmu_aalMU,1);
figure('Position',[0 0 350 1000]), heatmap([natfmu_aal natfmu_aalMU], 'YDisplayLabels',aal_label_reducMU)
colormap("jet_omega_mod")
