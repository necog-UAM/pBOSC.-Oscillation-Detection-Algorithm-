p.scripts = '\\150.244.69.206\nas\OMEGA\Enrique\OMEGA-Repo\Scripts';
addpath(genpath(p.scripts))

p.pBOSC = fileparts(p.scripts);
p.data = fullfile(p.pBOSC, 'OMEGA');
p.files = fullfile(p.pBOSC, 'Files');
p.raw = fullfile(p.pBOSC, 'Raw');

% Toolbox paths
addpath('Z:\Toolbox\fieldtrip-20230118') % Path to Fieldtrip
addpath('Z:\Toolbox\NECOG')
ft_defaults
ft_warning off

srate = 128;

% Participant and session list
[Nsub, subs, sess] = Omega_subs();
[source] = Omega_source(1925, 1.5);

% Frex parameters 
frex = exp(0.6:0.1:3.7); % 1.8 Hz to 40 Hz
Nfrex = length(frex);


Ntp = 33280 % CHANGE from OMEGA_pBOSC
Nvoxin = length(find(source.inside));
total_osc_time = zeros(Nsub,Nvoxin); % in all freqs
freq_osc_perc = zeros(Nsub,Nvoxin, length(frex)); % for each freq
freqbands_osc_perc = zeros(Nsub,Nvoxin, 5); % for each band

%%
tic
for s=1:Nsub
    fprintf('\nSubject:  %d / %d\n\n', s, Nsub)
    cd(['Z:\OMEGA\OMEGA_data\' 'sub-' subs{s} '\ses-' sess{s} '\conepis'])
    episname = 'episorig_v';
    for nepis = 1:Nvoxin 
        load([episname num2str(nepis)])
        Epis{nepis} = epis;
    end
    epis=Epis; %<--- rewrite better

%% Transform  epis from cells into a matrix

episodes = int8(zeros(Nvoxin,Nfrex,Ntp));
for v=1:Nvoxin
    epispnts = nan(length(epis{v}), length(frex));
    for ep=1:length(epis{v})
            fm = epis{v}(ep).freq;
            fbin = dsearchn(frex',fm);
            tpts = epis{v}(ep).timeps;
            episodes(v,fbin,tpts) = 1;
    end
end


%% Save individual cycles
episodes3 = zeros(size(episodes));
subepis_cycles = cell(size(episodes,1), 1);
for vx = 1:size(episodes,1)
    subepis_cycles{vx} = cell(size(episodes, 2), 1);
    for fx = 1:size(episodes,2)
            vector = squeeze(episodes(vx,fx, :))';
            % Pad edges
            pad = [0, vector, 0];
            starts = find(diff(pad) == 1);
            ends   = find(diff(pad) == -1) - 1;

            subepis_cycles{vx}{fx} = ends - starts + 1;
            subepis_cycles{vx}{fx} = (subepis_cycles{vx}{fx} ./ srate) * frex(fx);

            idx3cyc = subepis_cycles{vx}{fx} >= 3; % Index of episodes with at least 3 cycles
            tpts3cyc = arrayfun(@(s, e) s:1:e, starts(idx3cyc), ends(idx3cyc), 'UniformOutput', false); % Guardar los timepoints de >3 cycles
            for ep3=1:length(tpts3cyc)
                episodes3(vx,fx,tpts3cyc{ep3}) = 1;
            end
    end
end


save(['Z:\OMEGA\Enrique\Results\Cycles\nsmthcycles_s' num2str(s)], "subepis_cycles");
episodes = episodes3; % <-------- 3 cyles -----------

        % Total osc time is the proportion of time each voxel oscillates (any freq)
        tmp = squeeze(sum(episodes,2)>0);
        tmp(tmp>0) = 1;
        tmp = sum(tmp>0,2);
        total_osc_perc(s,:) = tmp ./ size(episodes,3) * 100;

        % Freq osc time is the proportion of time each voxel oscillates at each freq
        freq_osc_perc(s,:,:) = sum(episodes,3) ./  size(episodes,3) * 100; % size(episodes,3
        freq_osc_time_secs(s,:,:) = sum(episodes,3) ./ srate;
        freq_osc_time_cycles(s,:,:) = squeeze(freq_osc_time_secs(s,:,:)) .* frex;

        %      %% Freq bands osc time is the proportion of time each voxel oscillates at each band freq
        % for fb=1:5
        %     temp = squeeze(sum(squeeze(episodes(:,fb1(fb):fb2(fb),:)),2));
        %     temp(temp>0) = 1;
        %     freqbands_osc_perc(s,:,fb) = sum(temp,2) ./  size(episodes,3) * 100;%size(episodes,3);
        %     freqbands_osc_secs(s,:,fb) = sum(temp,2) ./ srate;%size(episodes,3);
        % end
end
t = toc;

cd('Z:\OMEGA\Enrique\Results')
save oscillatory_results_orig3cyc freq_osc_time_cycles freq_osc_time_secs freq_osc_perc total_osc_perc % freqbands_osc_secs freqbands_osc_perc
% load oscillatory_results_orig3cyc.mat % >3 cycles or >1 cycles

%% Percent (non)-oscillatory time

% VOXEL OSC TIME
        mean_total_osc_perc = squeeze(median(total_osc_perc,1))';
        pBOSC_nii( mean_total_osc_perc, 'total_osc_time')
        pBOSC_nii( 100-mean_total_osc_perc, 'total_nonosc_time')

% VOXEL-FREQ OSC TIME
        mean_freq_osc_perc = squeeze(mean(freq_osc_time_cycles,1));
        for fx =1:Nfrex
            tmpfreqperc = mean_freq_osc_perc(:,fx);
            pBOSC_nii(tmpfreqperc, [num2str(frex(fx)) '_orig'])
        end

% VOXEL-FREQ OSC CYC
mean_freq_osc_cyc = squeeze(mean(numbcycles,1, 'omitmissing'));
figure, bar(1:length(frex),mean(mean_freq_osc_cyc))
xticks(1:length(frex)), xticklabels(round(frex,1))

% VOXEL-FREQBAND OSC TIME


        mean_freqbands_osc_perc = squeeze(mean(freqbands_osc_perc,1));
        for fb = 1:length(fb_lab)
            dat = mean_freqbands_osc_perc(:,fb);
            Omega_nii(dat, ['mean_time_fband_' fb_lab{fb}])
        end
        mean_freqbands_osc_secs = squeeze(mean(freqbands_osc_secs,1));
        for fb = 1:length(fb_lab)
            dat = mean_freqbands_osc_secs(:,fb);
            Omega_nii(dat, ['new_mean_secs_fband_' fb_lab{fb}])
        end

 % Figures

 % NON-OSC TIME %
colcfg.colim   = [0 100];  
colcfg.colmap   = 'hot_omega_mod';
Omega_sourcefig(mean_total_osc_perc, 'osc', colcfg)

% Freqbands
fblims = [5 25; 5 25; 5 35; 5 25; 5 25];
fblabel = {'delta' 'theta' 'alpha' 'beta-low' 'beta-high'};

for fb = 1:5
colcfg.colim   = fblims(fb,:);  
colcfg.colmap   = 'hot_omega_mod';
Omega_sourcefig(mean_freqbands_osc_perc(:,fb), ['osctime_' fblabel{fb}], colcfg)
end

% All freqs

for fx = 1:size(mean_freq_osc_perc,2)
    tmpfreqperc = mean_freq_osc_perc(:,fx);
    if sum(tmpfreqperc) == 0
        continue
    else
        colcfg.colim   = [0 1];
        colcfg.colmap   = hot_omega_mod;
        Omega_sourcefig(tmpfreqperc, [num2str(fx) '_new_freq_osc_perc'], colcfg)
    end
end



%% Roi images
% ROIs
load(['Z:\OMEGA\OMEGA-NaturalFrequencies-main\mat_files\aal_voxel_label_10mm.mat'])
Nroi = length(aal_label_reduc);

% Roi index
for roi = 1:Nroi
    roiinds{roi} = find(label_inside_aal_reduc==roi);
end



% Roi aissign
for roi = 1:40
roifig = zeros(size(numbcycles,2),1);

    inds = find(label_inside_aal_reduc==roi);
    voxs = voxel_inside_aal(inds);
    roifig(voxs) = 1;
    Omega_nii(roifig, [num2str(roi) ' ' aal_label_reduc{roi}])

end

%% Oscillatory time
for rr = 1:Nroi
roimeans = squeeze(mean(mean_total_osc_perc(roiinds{rr}), 'omitmissing')); % select ROI to take mean cycs x frq
nosc = 100 - roimeans;
h = piechart([roimeans nosc]), title([aal_label_reduc{rr}]), 
saveas(gcf,[ num2str(rr) ' ' aal_label_reduc{rr}], 'png')
rr
close
end

%% Mean and sd duration of cycles
meancycles = NaN(Nsub,Nvoxin,Nfrex);
sdcycles = NaN(Nsub,Nvoxin,Nfrex);
epis_cycles = cell(Nsub,1);

for s = 1:Nsub
    s
    load(['Z:\OMEGA\Enrique\Results\Cycles\cycles_s' num2str(s)]); % original
    epis_cycles3cyc = subepis_cycles;
    for vx = 1:size(epis_cycles3cyc,1)
        meancycles(s,vx,:) = cellfun(@(x) mean(x(x>=3), 'omitmissing'), epis_cycles3cyc{vx});
        sdcycles(s,vx,:) = cellfun(@(x) std(x(x>=3), 'omitmissing'), epis_cycles3cyc{vx});
    end
end
save average_cycles meancycles sdcycles

gavgcycles = squeeze(mean(meancycles,1));
gsdcycles = squeeze(mean(sdcycles,1));
addpath(genpath('Z:\OMEGA\OMEGA-NaturalFrequencies-main\'))

% Figure
clmap = colormap(jet_omega_mod);
vx = 450;
avg = gavgcycles(vx,:);
sd = gsdcycles(vx,:);
groups = 1:32;              % group indices
% Define a colormap (32 distinct colors)
colors = jet_omega_mod;  % you can try 'parula(32)', 'jet(32)', etc.
colors = colors(2:2:end,:);

figure
hold on
for i = 1:32
    errorbar(avg(i), groups(i), sd(i), 'horizontal', 'o', ...
        'Color', colors(i,:), ...
        'MarkerFaceColor', colors(i,:), ...
        'CapSize', 6, 'LineWidth', 1.3);
end
set(gca,'YTick',groups), title([vx])
yticklabels(frex), xlim([0 20])

%% One general plot
oneavg = mean(gavgcycles, 'omitmissing');
onesd = mean(gsdcycles, 'omitmissing');
groups = 1:32;              % group indices

% Define a colormap (32 distinct colors)
colors = jet_omega_mod;  % you can try 'parula(32)', 'jet(32)', etc.
colors = colors(2:2:end,:);

figure
hold on

for i = 1:32
    errorbar(oneavg(i), groups(i), onesd(i), 'horizontal', 'o', ...
        'Color', colors(i,:), ...
        'MarkerFaceColor', colors(i,:), ...
        'CapSize', 6, 'LineWidth', 1.3);
end
set(gca,'YTick',groups)
yticklabels(frex)

%% El más característico de cada frecuencia
%
fblims = [4 8 14 22 33];
fblimsidx = dsearchn(frex', fblims');
fblabel = {'delta' 'theta' 'alpha' 'beta-low' 'beta-high'};
%
avg_epis_cycles_norm = normalize(avg_epis_cycles,1, "zscore");
idx = [];
for fx = 1:size(avg_epis_cycles_norm,2)
    [~,idx(fx)] = max(avg_epis_cycles_norm(:,fx));
end

clmap = colormap(jet_omega_mod);
groups = 1:32;
colors = jet_omega_mod;  % you can try 'parula(32)', 'jet(32)', etc.
colors = colors(2:2:end,:);
for vx = 2:31
    figure,
    avg = avg_epis_cycles(idx(vx),:);
    sd = sd_epis_cycles(idx(vx),:);
    hold on
    for i = 1:32
        errorbar(avg(i), groups(i), sd(i), 'horizontal', 'o', ...
            'Color', colors(i,:), ...
            'MarkerFaceColor', colors(i,:), ...
            'CapSize', 6, 'LineWidth', 1.3);
    end
    set(gca,'YTick',groups)
    yticklabels(frex)
    title([frex(vx) idx(vx)])
end

fblims = [4 8 14 22 33];
fblimsidx = dsearchn(frex', fblims');
fblabel = {'delta' 'theta' 'alpha' 'beta-low' 'beta-high'};
it=1;
for fb = 1:size(fblabel,2)
    dat = mean(avg_epis_cycles_norm(:,it:fblimsidx(fb)),2,'omitmissing');
    pBOSC_nii(dat,[fblabel{fb} '_cycles']);
    it = fblimsidx(fb);
end
