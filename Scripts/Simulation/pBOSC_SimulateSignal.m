%% Step 1. Simulate a dipole at the source level and project it to sensors. Aperiodic and gaussian noise are also added.

function simsignal = pBOSC_SimulateSignal(cfg)

%--------------------------------------------------------%

% cfg.signaldur = duration in seconds of the simulated signal [Default 60
% seconds].
% cfg.fsample = sample rate of the simulated signal [Default 512 Hz].
% cfg.dip = voxel where oscillations are generated.
% cfg.apgens = voxels where aperiodic component is generated.
% cfg.simcycles = cycles of each of the signals
% cfg.simfreq = frequency of each of the signals
% cfg.simsnr = signal to noise ratio of each of the signals
% cfg.figflag = generate figures (1) or not (0).
%--------------------------------------------------------%

%% Defaults
% Default duration of the signal is 60 seconds
if ~isfield(cfg, 'signaldur')
    signaldur = 60; % in seconds
else
    signaldur = cfg.signaldur;
end
% Default sample rate is 512
if ~isfield(cfg,'fsample')
    fsample = 512;
else
    fsample = cfg.fsample;
end
% Generate figures
if ~isfield(cfg, 'figflag')
    figflag = 0;
else
    figflag = cfg.figflag;
end


% Load source models
[source] = Omega_source(3423, 1.5);

dip = cfg.dip;
apgens = cfg.apgens;
simfreq = cfg.simfreq;
simcycles = cfg.simcycles;
simsnr = cfg.simsnr;

voxin = source.inside;
Nvox = length(voxin);

%% 1.1 Simulate signal

time  = 1/fsample:1/fsample:signaldur;
pnts  = length(time);

brainvox = zeros(Nvox, pnts); % voxels x time

% A. Aperiodic Generators
for gg = 1:length(apgens)
    % A. 1/f noise
    aperiodic = f_alpha_gaussian(pnts,1,1.2)';
    brainvox(apgens(gg),:) = aperiodic;
end

    % Figure: Aperiodic data
    if figflag 
        pinkcol = [100,1,100]/255;
        inside = find(source.inverse.inside);
        figure('Position',[0 0 1000 800])
        plot(source.inverse.pos(inside,2), source.inverse.pos(inside,3),'ok', 'MarkerSize',15, 'LineWidth',2)
        set(gca, 'Ylim', [-12 14])
        set(gca, 'Xlim', [-14 10])
        axis off
        for plt = 1:length(apgens)
        hold on
            plot(source.inverse.pos(inside(apgens(plt)),2), source.inverse.pos(inside(apgens(plt)),3),'ok', 'MarkerFaceColor', pinkcol, 'MarkerSize',19)
            pinkcol = pinkcol + 35/255;
        end
        pinkcol = [100,1,100]/255;
        ylm = [ -max(max(abs(brainvox))) max(max(abs(brainvox)))];
        for plt = 1:length(apgens)
        figure('Position', [0 500 1000 200])
        plot(brainvox(apgens(plt),:), 'Color', pinkcol, 'LineWidth',1)
        ylim(ylm), axis off
        pinkcol = pinkcol + 35/255;
        end
    end


% Leadfield orientations
lfvx = zeros(length(source.forward.grid.label),Nvox);
for vx = 1:Nvox
    lfvx(:,vx) = source.forward.grid.leadfield{voxin(vx)}(:,randi(3));
end

% Projection
noisedata = lfvx*brainvox;

    % Figure: Aperiodic in one sensor
    if figflag
    figure('Position', [0 500 1000 200])
    plot(noisedata(randi(size(noisedata,1)),:),'m', 'LineWidth',1), axis off
    axis off
    end

% Add gaussian noise on sensors 
gausnoise = randn(size(noisedata,1), size(noisedata,2)) .* prctile(prctile(abs(noisedata),50),50);
noisedata = noisedata + gausnoise;
noisedata_norm = noisedata ./ rms(noisedata,2);

    % Figure: White noise in one sensor
    if figflag
    ylm = [-max(max(abs(gausnoise))) max(max(abs(gausnoise)))];
    [~, idx] = max(rms(gausnoise,2))
    figure('Position', [0 500 1000 200])
    plot(gausnoise(idx,:),'Color',[128;128;128]/255, 'LineWidth',1), axis off
    ylim(ylm), axis off
    end

% B. Oscillation generators

brainvox = zeros(Nvox, pnts);
    cfg             =[];
    cfg.cycles      = simcycles;
    cfg.freq        = simfreq;
    cfg.fsample     = fsample;
    cfg.time        = time;
    oscillation     = Sim_Oscillation(cfg);
    brainvox(dip,:) = oscillation;
    osc_tps = find(oscillation); % save oscillation timepoints

    % Figure: Simulated oscillation
    if figflag
    figure('Position', [0 500 1000 200])
    plot(time, oscillation), axis off
    end


% Oscillator into leadfield
% Projection
signaldata = lfvx*brainvox;

% Find the channel with strongest projection
[~,chidx] = max(rms(signaldata,2));

    % Figure: Projected oscillation in one sensor
    if figflag
    figure('Position', [0 500 1000 200])
    plot(time, signaldata(chidx,:)), axis off
    end

% Normalize signal
signrms = rms(signaldata(chidx,osc_tps));
signaldata_norm = simsnr.*signaldata ./ signrms;

simdata = noisedata_norm + signaldata_norm;

    % Figure:
    if figflag
    ylm = [-max(max(abs(simdata))) max(max(abs(simdata)))];
    figure('Position', [0 500 1000 200])
    plot(time, simdata(chidx,:), 'k'), axis off
    end

% Store data in a struct
simsignal = [];
simsignal.cfg = [];
simsignal.fsample = fsample;
simsignal.trial{1} = simdata;
simsignal.time{1} = time;
simsignal.label = source.forward.grad.label;
simsignal.label = ft_channelselection('M*', simsignal); % select magnmtrs
simsignal.sim.snr = simsnr;
simsignal.sim.cycles = simcycles;
simsignal.sim.freq = simfreq;
simsignal.sim.dipole = dip;
simsignal.sim.numapgens = length(apgens);
simsignal.sim.apgens = apgens;
simsignal.sim.osctimepoints = osc_tps;
simsignal.sim.fsample = fsample;

end