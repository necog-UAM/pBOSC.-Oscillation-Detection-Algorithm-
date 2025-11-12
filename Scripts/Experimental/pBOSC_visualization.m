%% 0. Init
[Nsub, subs, sess] = Omega_subs();

dpath = 'Z:\OMEGA\OMEGA_data\';
rawpath = 'Z:\OMEGA\OMEGA_raw\';       
addpath('Z:\OMEGA\Enrique')
addpath(genpath('Z:\Toolbox'))
addpath('C:\Users\Usuario\Desktop\WorkingDir\OMEGA_Repo\Scripts\BaseScripts')

s = 2;
sub = subs{s};
ses = sess{s};
       
%% Compare noise thresholds
cd([dpath '\sub-' sub '\ses-' ses])
load thshld
load thshld_noise

for vx = 1:size(thshld,1)
plot(thshld(vx,:))
hold on, plot(thshld_noise(vx,:))
pause(0.1)
clf
end

%% Load data and episodes

load datasource_allvox.mat
load aperiodic_data.mat

% Correction
dataclean.trial{1} = dataclean.trial{1} ./ rms(aperiodic,2);
noise_corr = aperiodic ./ rms(aperiodic,2);

cfg = [];
cfg.begsample = 1;
cfg.endsample = length(aperiodic);
dataclean = ft_redefinetrial(cfg,dataclean)
    
% IV.b Analyze Whole signal
dataclean_ap = dataclean;
dataclean_ap.trial{1} = noise_corr;
%
% V. Downsample both signals
cfg = [];
cfg.resamplefs = 128;
[dataclean] = ft_resampledata(cfg, dataclean);
[dataclean_ap] = ft_resampledata(cfg, dataclean_ap);
fsample = dataclean.fsample;

   %% Step 5.2.Frequency analysis
  
   [frex, Nfrex, fb_lab,fb1,fb2] = Omega_frex();


   % I. Original signal
   warning('off')
   Nvox = length(dataclean.label);
   Ntp = length(dataclean.time{1});
   powspctm = single(zeros(Nvox, Nfrex));
   powspctm = repmat(powspctm, [1 1 Ntp]);

   for f0 = 1:Nfrex
       cfg            = [];
       cfg.method     = 'mtmconvol';
       cfg.taper      = 'hanning';
       cfg.output     = 'pow';
       cfg.foi        = frex(f0);
       cfg.toi        = dataclean.time{1};
       cfg.t_ftimwin  = 5./frex(f0);
       cfg.pad        = 'nextpow2';
       cfg.keeptrials = 'yes';
       freq           = ft_freqanalysis(cfg, dataclean);

       powspctm(:,f0,:)= single(squeeze(freq.powspctrm));
   end

   % II. Aperiodic signal
   powspctmap = single(zeros(Nvox, Nfrex));
   powspctmap = repmat(powspctmap, [1 1 Ntp]);

   for f0 = 1:Nfrex
       cfg            = [];
       cfg.method     = 'mtmconvol';
       cfg.taper      = 'hanning';
       cfg.output     = 'pow';
       cfg.foi        = frex(f0);
       cfg.toi        = dataclean_ap.time{1};
       cfg.t_ftimwin  = 5./frex(f0);
       cfg.pad        = 'nextpow2';
       cfg.keeptrials = 'yes';
       freqap           = ft_freqanalysis(cfg, dataclean_ap);

       powspctmap(:,f0,:)= single(squeeze(freqap.powspctrm));
   end
   

       % % Plot histograms of threshold
    %     vxplot = voxin(1357); fx = dsearchn(frex',10);
    %     figure, h = histogram(squeeze(powspctm(vxplot,fx,:)),100),
    %     hold on, histogram(squeeze(powspctmap(vxplot,fx,:)),'BinWidth',h.BinWidth),
    %     hold on, xline(thshld(vxplot,fx), 'r--', 'LineWidth',3)
    %     title([num2str(frex(fx)) 'Hz, ' '| Voxel ' num2str(vxplot)]), legend({'Signal','Aperiodic'})
    %     text(thshld(fx)+thshld(fx)/2, mean(ylim), num2str(thshld(fx)), 'FontSize',16)
    %     %
    
        vx = 665; fx = dsearchn(frex',3);
            figure
            subplot(211), h = histogram(squeeze(powspctm(vx,fx,:)),100),    
            hold on, histogram(squeeze(powspctmap(vx,fx,:)),'BinWidth',h.BinWidth),
            xl = xlim();
            % hold on, xline(thshld95(fx), 'r--', 'LineWidth',3)
            % hold on, xline(thshld99(fx), 'r--', 'LineWidth',3)
            % % hold on, xline(vxfxths(vx,fx), 'k--', 'LineWidth',3)
            % hold on, xline(vxfxths2(vx,fx), 'm--', 'LineWidth',3)
            % hold on, xline(vxfxths3(vx,fx), 'g--', 'LineWidth',3)
            hold on, xline(thshld(vx,fx), 'g--', 'LineWidth',3)
            xline(thshld_noise(vx,fx) , 'k--', 'LineWidth',3)

            subplot(212),  hold on
            % plot(cntrs(fx,2:end), snratio1(fx,:),'.')
            % plot(cntrs(fx,:), 1-snratio2(fx,:),'m.')
            % plot(cntrs(fx, idx(fx)), 1-snratio2(fx,idx(fx)),'ms','LineWidth',2)
            % plot(cntrs2(fx,:), 1-snratio3(fx,:),'g*')
            % plot(cntrs2(fx, idx2(fx)), 1-snratio3(fx,idx2(fx)),'gs','LineWidth',2)
            plot(cntrs(fx,:), snratio(fx,:),'b.')
            plot(cntrs(fx, idx(fx)), snratio(fx,idx(fx)),'rs','LineWidth',2)
            yline(.95, '--', 'LineWidth',2)
            xlim(xl)

            subplot(211), xlim([.1 .4])
            subplot(212), xlim([.1 .4])


 %% -- Episodes

[voxel_inside, Nvox, voxin, Nvoxin] = Omega_VoxInside();

    frex = [0 frex]
    Nsub = 1;

    epis_number_all = zeros(Nsub, Nvoxin ,length(frex));
    epis_power_all = zeros(Nsub, Nvoxin ,length(frex));
    epis_cycles_all = zeros(Nsub, Nvoxin ,length(frex));

for s=1:Nsub
cd([dpath '\sub-' subs{s} '\ses-' sess{s}])
load epis_prctile95noise_1cyc.mat
episodes = epis;
clear epis
   
epis_number = zeros(length(episodes),3000);
epis_power = zeros(length(episodes),length(frex));
epis_cycles = zeros(length(episodes),length(frex));

for ii = 1:length(episodes)

    if isempty(episodes{ii})
        episodes{ii}.freq = 0;
        episodes{ii}.dur_cyc = 0;
        episodes{ii}.power = 0;
    end

epis_number(ii,:) = [[episodes{ii}.freq] zeros(1,size(epis_number,2)- length([episodes{ii}.freq]))];

[val,~,ix] = unique(frex(dsearchn(frex',epis_number(ii,:)')'));
episcounts = accumarray(ix,1)';
[~,idx] = ismember(val,frex);

tmpv = num2cell(frex(idx(ix(1:length([episodes{ii}.freq])))));
[episodes{ii}.cntrs] = tmpv{:};

for lp = 1:length(frex)
    idx2 = [episodes{ii}.cntrs] == frex(lp);
    epis_power(ii,lp) = sum([episodes{ii}(idx2).power]); % power
    epis_cycles(ii,lp) = sum([episodes{ii}(idx2).dur_cyc]); % cycles
end

epis_number_all(s,ii,idx) = episcounts; % counts
epis_power_all(s,:,:) = epis_power; 
epis_cycles_all(s,:,:) = epis_cycles; 

end
s
end


% Plots
epis_number_norm = normalize(epis_number_all,2,"zscore", "robust");
epis_number_mean = squeeze(mean(epis_number_norm(:,:,2:end),1,'omitmissing'));

epis_power_norm = normalize(epis_power_all,2,"zscore", "robust");
epis_power_mean = squeeze(mean(epis_power_norm(:,:,2:end),1,'omitmissing'));

epis_cycles_norm = normalize(epis_cycles_all,2,"zscore", "robust");
epis_cycles_mean = squeeze(mean(epis_cycles_norm(:,:,2:end),1,'omitmissing'));

voxlab = {'Hippcmps', 'Mtr', 'Frntl', 'Prcneus', 'ACC' 'Out' 'Center', 'Partl'};
voxls = [298 1829 759 1357 1468 1794 660 1384];
figure(1),
for vx = 1:length(voxls);
    subplot(2,4,vx)
    bar(epis_number_mean(voxls(vx),:),'k'), xticks(2:length(frex)), xticklabels(frex(2:end))
    title(voxlab{vx}), sgtitle('Number of episodes')
end

figure(2),
for vx = 1:length(voxls);
    subplot(2,4,vx)
    bar(epis_power_mean(voxls(vx),:),'b'), xticks(2:length(frex)), xticklabels(frex(2:end))
    title(voxlab{vx}), sgtitle('Power')
end

figure(3),
for vx = 1:length(voxls);
    subplot(2,4,vx)
    bar(epis_cycles_mean(voxls(vx),:),'r'), xticks(2:length(frex)), xticklabels(frex(2:end))
    title(voxlab{vx}), sgtitle('Cycles')
end

