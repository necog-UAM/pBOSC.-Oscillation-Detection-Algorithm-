%% Evaluate results

    clear all
    
    dpath = 'Z:\OMEGA\OMEGA_data\';
    addpath('Z:\Toolbox\fieldtrip-20230118')
    addpath(genpath('C:\Users\Usuario\Desktop\WorkingDir\OMEGA_Repo\Scripts'))
    
    %% Step 0. Subjects
    [Nsub, subs, sess] = Omega_subs();
    
    %% Step 0. Sources
    [voxel_inside, Nvox, voxin, Nvoxin] = Omega_VoxInside();
    
    %% Step 0.Frequencies
    [frex, Nfrex, fb_lab,fb1,fb2] = Omega_frex();
    
    %% Step 0.Define neighbouring voxels
    [dim,xx,yy,zz,connmat, dtempl] = Omega_neighbors(1.5);
    
    %% Explore epis
    frex = [0 frex]
    Nsub = 20;

    epis_number_all = nan(Nsub, Nvoxin ,length(frex));
    epis_power_all = nan(Nsub, Nvoxin ,length(frex));
    epis_cycles_all = nan(Nsub, Nvoxin ,length(frex));
    thshld_noise = nan(Nsub, Nvoxin, length(frex)-1);

    for s=2:Nsub
        cd([dpath 'sub-' subs{s} '\ses-' sess{s}])
if exist('epis.mat')==2
        load epis.mat
        load thshld.mat

        epis_number = zeros(length(epis),5000);
        epis_power = zeros(length(epis),length(frex));
        epis_cycles = zeros(length(epis),length(frex));

        for ii = 1:length(epis)

            if isempty(epis{ii})
                epis{ii}.freq = 0;
                epis{ii}.dur_cyc = 0;
                epis{ii}.power = 0;
            end

            epis_number(ii,:) = [[epis{ii}.freq] zeros(1,size(epis_number,2)- length([epis{ii}.freq]))];

            [val,~,ix] = unique(frex(dsearchn(frex',epis_number(ii,:)')'));
            episcounts = accumarray(ix,1)';
            [~,idx] = ismember(val,frex);

            tmpv = num2cell(frex(idx(ix(1:length([epis{ii}.freq])))));
            [epis{ii}.cntrs] = tmpv{:};

            for lp = 1:length(frex)
                idx2 = [epis{ii}.cntrs] == frex(lp);
                epis_power(ii,lp) = sum([epis{ii}(idx2).power]); % power
                epis_cycles(ii,lp) = sum([epis{ii}(idx2).dur_cyc]); % cycles
            end

            epis_number_all(s,ii,idx) = episcounts; % counts
            epis_power_all(s,:,:) = epis_power;
            epis_cycles_all(s,:,:) = epis_cycles;
            thshld_noise(s,:,:) = thshld(voxin,:);
        end
end
    end

 


% Plots
epis_number_norm = normalize(epis_number_all,2,"zscore", "robust");
epis_number_mean = squeeze(mean(epis_number_norm(:,:,2:end),1,'omitmissing'));

epis_power_norm = normalize(epis_power_all,2,"zscore", "robust");
epis_power_mean = squeeze(mean(epis_power_norm(:,:,2:end),1,'omitmissing'));

epis_cycles_norm = normalize(epis_cycles_all,2,"zscore", "robust");
epis_cycles_mean = squeeze(mean(epis_cycles_norm(:,:,2:end),1,'omitmissing'));

%% EPISODES PLOTS
voxlab = {'Hippcmps', 'Mtr', 'Frntl', 'Prcneus', 'CRBLM' 'Out' 'Center', 'Partl'};
voxls = [148 1739 558 1357 413 1794 1067 1703];
figure(1),
for vx = 1:length(voxls);
    subplot(2,4,vx)
    bar(epis_number_mean(voxls(vx),:),'b', 'FaceAlpha',0.5), xticks(2:length(frex)), xticklabels(frex(2:end))
    title(voxlab{vx}), sgtitle('Number of episodes ')
end

figure(2),
for vx = 1:length(voxls);
    subplot(2,4,vx)
    bar(epis_power_mean(voxls(vx),:),'b', 'FaceAlpha',0.5), xticks(2:length(frex)), xticklabels(frex(2:end))
    title(voxlab{vx}), sgtitle('Power')
end

figure(3),
for vx = 1:length(voxls);
    subplot(2,4,vx)
    bar(epis_cycles_mean(voxls(vx),:),'b', 'FaceAlpha',0.5), xticks(2:length(frex)), xticklabels(frex(2:end))
    title(voxlab{vx}), sgtitle('Cycles')
end

%% Thresholds plots
figure
for vx = 1:length(voxls)
    subplot(3,3,vx)
    plot(mean(squeeze(thshld(:,voxls(vx),:)))), hold on
    plot(mean(squeeze(thshld_signnoise(:,voxls(vx),:))))
    title(voxlab{vx}), sgtitle('Thresholds')
end




%% All subjects analysis

Nsub = 20;
dur_cyc = zeros(Nsub,Nvoxin,Nfrex);
dur_sec = zeros(Nsub,Nvoxin,Nfrex);
num_epis = zeros(Nsub,Nvoxin,Nfrex);
mean_power = zeros(Nsub,Nvoxin,Nfrex);

ct=1;
for s=2:Nsub
    ct
    cd([dpath '\sub-' subs{s} '\ses-' sess{s}])
    if exist('epis.mat')==2
    load(['epis.mat'])

    for v=1:Nvoxin
        for ep=1:length(epis{v})
            % if epis{v}(ep).dur_cyc >=2
            fm = epis{v}(ep).freq;
            fbin = findbin(frex,fm);
            dur_cyc(ct,v,fbin)=epis{v}(ep).dur_cyc + dur_cyc(ct,v,fbin);
            dur_sec(ct,v,fbin)=epis{v}(ep).dur_sec + dur_sec(ct,v,fbin);
            num_epis(ct,v,fbin)=num_epis(ct,v,fbin)+1;
            mean_power(ct,v,fbin)=epis{v}(ep).power + mean_power(ct,v,fbin);
            % end
        end
    end
    mean_power2(ct,:,:) = mean_power(ct,:,:)./num_epis(ct,:,:);
    dur_cyc2(ct,:,:) = dur_cyc(ct,:,:)./num_epis(ct,:,:);
    ct = ct+1;
    end
end

% mean_power(isnan(mean_power))=0;
% dur_cyc(isnan(dur_cyc))=0;

dur_cycm = squeeze(nanmean(dur_cyc,1));
dur_secm = squeeze(nanmean(dur_sec,1));
num_epism = squeeze(nanmean(num_epis,1));
mean_powerm = squeeze(nanmean(mean_power,1));
mean_powerm2 = squeeze(nanmean(mean_power2,1));
dur_cycm2 = squeeze(nanmean(dur_cyc2,1));



% Load
cd([dpath 'sub-' subs{1} '\ses-' sess{1}])
load('source_inverse_10mm.mat')
voxel_inside = find(source.inside==1);
load('source_forward_10mm_allvox.mat')

% number of cycles/freq band
for f0 = 1:length(fb1)
    source2 = source;
    source2.avg.pow(voxel_inside) = nanmean(mean_powerm2(:,fb1(f0):fb2(f0)),2);
    source2.avg.mom=cell(length(source2.avg.noise),1);
    source2.time=1;

    cfg = [];
    cfg.parameter  = 'avg.pow';
    cfg.downsample = 2;
    cfg.interpmethod  =  'linear';
    source_interp  = ft_sourceinterpolate (cfg, source2, source_forward.mri);

    cd('C:\Users\Usuario\Desktop\WorkingDir\OMEGA_Repo\Results')% cd([dpath '\fBOSC'])
    cfg = [];
    cfg.filetype  = 'nifti';
    cfg.parameter = 'pow';
    cfg.filename  = [fb_lab{f0} 'epis3cyc-pow'];
    ft_sourcewrite(cfg, source_interp)
end

% Natfrex
all_nat = normalize(mean_power,2,"zscore");
all_nat = squeeze(nanmedian(all_nat,1));
% avg_rzscore = squeeze(nanmean(all_nat,1));
% avg_rzscore = movmean(avg_rzscore,5,1);     % +/-2 timepoints

nf=NaN(Nvoxin,1);
for v = 1:Nvoxin
    [pks,locs] = findpeaks(all_nat(v,:));
    try
        pks2 = pks(~isinf(pks));
        locs = locs(~isinf(pks));
        temp=frex(locs(find(pks2==max(pks2))));
        nf(v)=temp(1);
    end
end

% all


source2 = source;
source2.avg.pow(voxel_inside) = nf;
source2.avg.mom=cell(length(source2.avg.pow),1);
source2.time=1;

cfg = [];
cfg.parameter  = 'avg.pow';
cfg.downsample = 2;
cfg.interpmethod  =  'nearest';
source_interp  = ft_sourceinterpolate (cfg, source2, source_forward.mri);

cd(['C:\Users\Usuario\Desktop\WorkingDir\OMEGA_Repo\Results'])
cfg = [];
cfg.filetype  = 'nifti';
cfg.parameter = 'pow';
cfg.filename  = ['Natfrex_only-noise'];
ft_sourcewrite(cfg, source_interp)




%%%%%

Nvox = length(voxel_inside);

% dominant frequency
fdom = zeros(Nvox,1);
for v = 1:Nvox
    try
        fdom(v) = frex(find(dur_secm(v,1:58)==max(dur_secm(v,1:58))));
    end
end
dur_cycz = zscore(dur_secm,1);
fnat = zeros(Nvox,1);
for v = 1:Nvox
    try
        fnat(v) = frex(find(dur_cycz(v,1:58)==max(dur_cycz(v,1:58))));
    end
end

source2 = source;
source2.avg.pow(voxel_inside) = drms;
source2.avg.mom=cell(length(source2.avg.noise),1);
source2.time=1;

cfg = [];
cfg.parameter  = 'avg.pow';
cfg.downsample = 2;
cfg.interpmethod  =  'linear';
source_interp  = ft_sourceinterpolate (cfg, source2, source_forward.mri);

cd([naspath '\fBOSC'])
cfg = [];
cfg.filetype  = 'nifti';
cfg.parameter = 'pow';
cfg.filename  = ['sub1_rms_uncorr'];
ft_sourcewrite(cfg, source_interp)


source2 = source;
source2.avg.pow(voxel_inside) = log(fdom);
source2.avg.mom=cell(length(source2.avg.noise),1);
source2.time=1;

cfg = [];
cfg.parameter  = 'avg.pow';
cfg.downsample = 2;
cfg.interpmethod  =  'nearest';
source_interp  = ft_sourceinterpolate (cfg, source2, source_forward.mri);


figure
cfg               = [];
cfg.figure        = 'gca';
cfg.method        = 'surface';
cfg.funparameter  = 'pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [0.7 3.4];      % exponential scale   2-30 Hz
cfg.funcolormap   = 'jet_omega_mod';
cfg.projmethod    = 'nearest';
cfg.opacity       = 0.8;
cfg.camlight      = 'no';
cfg.colorbar      = 'yes';
cfg.surffile     = 'surface_pial_left.mat';
cfg.surfinflated  = 'surface_inflated_left_caret_white.mat';
subplot(2,2,1), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('left')
subplot(2,2,3), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('left')

cfg.surffile     = 'surface_pial_right.mat';
cfg.surfinflated  = 'surface_inflated_right_caret_white.mat';
subplot(2,2,2), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('right')
subplot(2,2,4), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('right')



% Dominant frequency in ROIs

cd('I:\Mis documentos\PROYECTOS\OMEGA\Kmeans\Nk25_10mm_tot\')
load centroid_labels

Nroi = max(centroid_labels);

f1 = 0;
f2 = 4.6;
d  = 0.2;
foi1 = exp(f1:d:f2);
foi2 = exp(f1:d/10:f2);
for i = 1:length(foi2)
    c0 = foi2(i);
    rectp(i,:) = rectangularPulse(c0-c0/4 ,c0+c0/4 ,foi2);      % to isolate individual peaks later
end
p = parpool(4);                  % parallel pool
p.IdleTimeout = inf;

natf_roi = [];
for i=1:Nroi
    ff = fdom(centroid_labels==i);

    [counts,~] = histcounts(ff,foi1);
    centers = exp(log(foi1) + d/2);
    counts  = interp1(centers(1:end-1),counts,foi2,'pchip');
    centers = foi2;

    [pks,locs] = findpeaks(counts);
    [~,ii] = sort(pks,'descend');
    if length(ii) >= 3
        c0s = centers(locs(ii(1:3)));         % fit a maximum of 3 peaks to speed up
    else
        c0s = centers(locs);
    end

    % gaussian fit of all candidate peak frequencies
    [A,mu,sigma,rsq] = gausfitc0(c0s,counts,centers,rectp,foi2);
    natf_roi(i) = mu(find(A==max(A)));
end

natf_vox = [];
for i=1:Nroi
    natf_vox(centroid_labels==i) = natf_roi(i);
end


source2 = source;
source2.avg.pow(voxel_inside) = log(natf_vox);
source2.avg.mom=cell(length(source2.avg.noise),1);
source2.time=1;

cfg = [];
cfg.parameter  = 'avg.pow';
cfg.downsample = 2;
cfg.interpmethod  =  'nearest';
source_interp  = ft_sourceinterpolate (cfg, source2, source_forward.mri);


figure
cfg               = [];
cfg.figure        = 'gca';
cfg.method        = 'surface';
cfg.funparameter  = 'pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [0.7 3.4];      % exponential scale   2-30 Hz
cfg.funcolormap   = 'jet_omega_mod';
cfg.projmethod    = 'nearest';
cfg.opacity       = 0.8;
cfg.camlight      = 'no';
cfg.colorbar      = 'yes';
cfg.surffile     = 'surface_pial_left.mat';
cfg.surfinflated  = 'surface_inflated_left_caret_white.mat';
subplot(2,2,1), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('left')
subplot(2,2,3), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('left')

cfg.surffile     = 'surface_pial_right.mat';
cfg.surfinflated  = 'surface_inflated_right_caret_white.mat';
subplot(2,2,2), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('right')
subplot(2,2,4), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('right')
