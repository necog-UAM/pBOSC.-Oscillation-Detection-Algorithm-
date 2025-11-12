% Omega Natural frequencies
addpath(genpath('Z:\OMEGA\Enrique\OMEGA-Repo\Scripts\'))
addpath Z:\Toolbox\fieldtrip-20230118
ft_defaults
rawpath = 'Z:\OMEGA\OMEGA_raw\'; 
dpath = 'Z:\OMEGA\OMEGA_data\'
addpath Z:\OMEGA\OMEGA-NaturalFrequencies-main\subfxs
ft_warning off
addpath('Z:\Toolbox\NECOG')

[frex, Nfrex, fb_lab,fb1,fb2] = Omega_frex();
[~, ~,~, ~,~,~,~,connmat] = Omega_voxin(1925, 1.5) ;

cd('Z:\OMEGA\Enrique\Results')
load natfreq_Nk25.mat
load oscillatory3_results_singlepis.mat

Nvoxin = size(freq_osc_time_cycles,2);
Nsub = size(freq_osc_time_cycles,1)
freq_osc_time_cycles = permute(freq_osc_time_cycles, [3, 2, 1]); %frex x vox x subs

%% Natural Frequencies of each participant
freqvoxs = zeros(size(freq_osc_time_cycles,2), size(freq_osc_time_cycles,3));
for s = 1:Nsub
    s
    single_frexvox = freq_osc_time_cycles(:,:,s);
    % Normalize across voxels
    single_frexvox(single_frexvox==0)=NaN;
    single_frexvox = normalize(single_frexvox,2, "zscore");
    for v = 1:Nvoxin
        vxneigh = connmat(vx,:)==1;
        freqvoxs_smth= single_frexvox(:,vxneigh);
        propinterp = [];
       for i=1:size(freqvoxs_smth,2)
            propinterp(:,i) = interp(freqvoxs_smth(:,i),10);
        end
        [h,p,ci,stats] = ttest(propinterp');

        [tmax,idmax] = max(stats.tstat);
        fnatsm(vx) = f2(idmax);
        fnatsm_T(vx) = tmax;
        fnatsm_p(vx) = p(idmax);

        tval = stats.tstat;
        tval(p>.05) = NaN;
        tval(tval < 0) = NaN;
        [pks,locs] = findpeaks(tval);

        fnat.pks{vx}=f2(locs);
        fnat.w{vx}=pks*100./sum(pks);

    end













%% Gaussian fit

% foi2 = exp(0.6:0.02:3.7);
f1 = 0;
f2 = 3.7;
d  = 0.2;

foi1 = exp(f1:d:f2);
foi2 = exp(f1:d/10:f2);

for i = 1:length(foi2)
    c0 = foi2(i);
    rectp(i,:) = rectangularPulse(c0-c0/4 ,c0+c0/4 ,foi2);      % to isolate individual peaks later
end

p = parpool(4);                  % parallel pool
p.IdleTimeout = inf;

Nvox = size(freqvoxs,1);
Nboot = 500;
natf.A       = NaN(Nvox,Nboot);
natf.mu      = NaN(Nvox,Nboot);
natf.sigma   = NaN(Nvox,Nboot);
natf.rsq     = NaN(Nvox,Nboot);
natf.randsub = NaN(Nboot,Nsub);

% bootstrap confidence interval; based on bootci -> ci = bootci(Nboot,{@bootfun,freqvoxs},'type','per')
rng('shuffle')
for b = 1:Nboot
    disp(['Bootstrapping ' num2str(b) '/' num2str(Nboot)])
    randsub = randi(Nsub,Nsub);
    randsub = randsub(1,:);  % sample with replacement Nsub
    for v = 1:Nvox
        fv = freqvoxs(v,randsub);
        [counts,~] = histcounts(fv(:),foi1);
        centers = exp(log(foi1) + d/2);
        figure(1), clf, plot(centers(1:end-1), counts,'--')        
       
        
        counts  = interp1(centers(1:end-1),counts,foi2,'pchip');
        centers = foi2;

        [pks,locs] = findpeaks(counts);
        [~,ii] = sort(pks,'descend');
        if length(ii) >= 3
            c0s = centers(locs(ii(1:3))); % fit a maximum of 3 peaks to speed up
        else
            c0s = centers(locs);
        end

        figure(1), hold on, plot(centers, counts)
        % gaussian fit of all candidate peak frequencies
        [A,mu,sigma,rsq] = gausfitc0_serial(c0s,counts,centers,rectp,foi2);

        i = find(rsq == max(rsq));% identify the one with the highest goodness of fit
        % i = find(A == max(A));% identify the one with the highest goodness of fit
        natf.A(v,b)       = A(i);
        natf.mu(v,b)      = mu(i);
        natf.sigma(v,b)   = sigma(i);
        natf.rsq(v,b)     = rsq(i);
        natf.randsub(b,:) = randsub;
        gausf = A(i) * exp(-(centers-mu(i)).^2 /(2*sigma(i).^2));
        plot(centers,squeeze(gausf),'k', 'LineWidth',2)
        set(gca,'XLim',[0 35],'YLim',[0 80])
        hold off
        pause
    end
end

natfmu      = median(natf.mu,2,"omitmissing");
natfstd = std(natf.mu,[],2, "omitmissing");
natfci(:,1) = prctile(natf.mu,2.5,2);          % 95% CI interval
natfci(:,2) = prctile(natf.mu,97.5,2);



cd('Z:\OMEGA\Enrique\Results')
% load natfreqgauss3
save natfreqgauss3 natf natfmu natfci

% Figures
Omega_nii(natfmu, 'nf_3_gaussfit')

colcfg = [];
colcfg.colim   = [0.7 3.4];  
colcfg.colmap   = 'jet_omega_mod';
Omega_sourcefig(log(natfmu), 'test',colcfg)

    % randsub = [1:16];
    % v=250;
    % fv = freqvoxs(v,randsub);
    % [counts,~] = histcounts(fv(:),frex);
    % centers = exp(log(frex) + d/2);
    % counts  = interp1(centers(1:end-1),counts,foi2,'pchip');
    % centers = foi2;
    % figure,plot(centers,counts)



%% Natural frequencies ROIs

% Load roi str
load(['Z:\OMEGA\OMEGA-NaturalFrequencies-main\mat_files\aal_voxel_label_10mm.mat'])
Nroi = length(aal_label_reduc);

% figure
for roi = 1:Nroi
    inds = find(label_inside_aal_reduc==roi);
    voxs = voxel_inside_aal(inds);
    muvoxs = natf.mu(voxs,:);
    
    Nk = 4;
    [idx,C,sumd,D] = kmeans(muvoxs,Nk,'Replicates',5,'MaxIter',200);      % k-means 4 clusters max
    nvoxs = [sum(idx==1) sum(idx==2) sum(idx==3) sum(idx==4) ];           % number of voxels in each cluster  
  Nk = Nk - sum(nvoxs==1);                                              % eliminate clusters with only 1 member and repeat clustering                                           
    [idx,C,sumd,D] = kmeans(muvoxs,Nk,'Replicates',5,'MaxIter',200);      % k-means adapting number of clusters 
    nvoxs = [sum(idx==1) sum(idx==2) sum(idx==3) sum(idx==4) ];           % number of voxels in each cluster                       
    dominclust = find(nvoxs == max(nvoxs));                               % dominant cluster
    centroidvox = find(D(:,dominclust)==min(D(:,dominclust)));            % voxel closer to centroid
    represvox(roi) = voxs(centroidvox(1));                                % representative voxel for this AAL region
     % hold on, histogram(C'),title(aal_label_reduc{roi}),pause                    % visualize pattern of natural frequencies for both clusters
end

natfmu_aal      = median(natf.mu(represvox,:),2);
natfci_aal(:,1) = prctile(natf.mu(represvox,:),2.5,2);       % 95% CI interval
natfci_aal(:,2) = prctile(natf.mu(represvox,:),97.5,2);

% load natfreq_aal
save natfreq_aal natfmu_aal natfci_aal represvox aal_label_reduc natf

figure, plot(natfci_aal)

% Roi aissign
roinf = zeros(1925,1);
for roi = 1:Nroi
    inds = find(label_inside_aal_reduc==roi);
    voxs = voxel_inside_aal(inds);
    roinf(voxs) = natfmu_aal(roi);
end

Omega_nii(roinf, 'Roi_NatfreqKmeans')

% Hist allvox per ROI
figure,
sgtitle('All vox')
for rr = 1:Nroi
    subplot(4,10,rr)
% boxplot(natf_aal(rr,:))
inds = find(label_inside_aal_reduc==rr);
voxs = voxel_inside_aal(inds);
histogram(natf.mu(voxs,:),frex)
xlim([0 35]), ylims(rr,:) = ylim();
title(aal_label_reduc{rr},'FontSize',8),
end

% Distribuciones
% Hist representative vox per ROI
figure,
sgtitle('Representative vox')
for rr = 1:Nroi
    subplot(4,10,rr)
% boxplot(natf_aal(rr,:))
histogram(natf.mu(represvox(rr),:),frex)
xlim([0 35])
title(aal_label_reduc{rr},'FontSize',8),
end



% Comparar
load Z:\OMEGA\Enrique\Results\natfreq_aalMU.mat
 fMU   = 0.55:0.05:4.6;
 foiMU = exp(fMU);
% Distribuciones de ALMU
% Hist allvox per ROI
figure,
sgtitle('All vox Almu')
for rr = 1:Nroi
    subplot(4,10,rr)
% boxplot(natf_aal(rr,:))
inds = find(label_inside_aal_reduc==rr);
voxs = voxel_inside_aal(inds);
histogram(natfMU.mu(voxs,:),frex)
xlim([0 35]), ylims(rr,:) = ylim();
title(aal_label_reduc{rr},'FontSize',8),
end

% Hist representative vox per ROI
figure,
for rr = 1:Nroi
    subplot(4,10,rr)
histogram(natfMU.mu(represvoxMU(rr),:), frex,'EdgeColor','none');
hold on
histogram(natf.mu(represvox(rr),:),frex,'EdgeColor','none');
xlim([0 35])
title(aal_label_reduc{rr},'FontSize',8),
end


% Violins
figure, 
for rr = 1:Nroi
    subplot(4,10,rr)
    violinplot(natfMU.mu(represvoxMU(rr),:), 'DensityDirection','negative')
    hold on
    violinplot(natf.mu(represvox(rr),:), 'DensityDirection','positive')
    title(aal_label_reduc{rr},'FontSize',8),
    ylim([0 30])
end

% Histograms

figure(1), 
for rr = 1:Nroi
    subplot(4,10,rr)
    histogram(natfMU.mu(represvoxMU(rr),:),log(frex))
    xlim([ 0.6 3.7])
    % hold on
    % histogram(log10(natf.mu(represvox(rr),:)),frex)
    % title(aal_label_reduc{rr},'FontSize',8),
end


% MriCron ALmu
% Roi aissign
roinfMU = zeros(1925,1);
for roi = 1:Nroi
    inds = find(label_inside_aal_reduc==roi);
    voxs = voxel_inside_aal(inds);
    roinfMU(voxs) = natfmu_aalMU(roi);
end

Omega_nii(roinfMU, 'Roi_NatfreqKmeans_ALMU')

natfmu_Almu      = median(natfMU.mu,2,"omitmissing");
colcfg = [];
colcfg.colim   = [0.7 3.4];  
colcfg.colmap   = 'jet_omega_mod';
Omega_sourcefig(log(natfmu_Almu), 'test',colcfg)


%% Correlacion
natfstd_Almu = std(natfMU.mu,[],2, "omitmissing");
natfci_Almu(:,1) = prctile(natfMU.mu,2.5,2);          % 95% CI interval
natfci_Almu(:,2) = prctile(natfMU.mu,97.5,2);

% Corr voxels
[rho,  pval] = corr(natfmu, natfmu_Almu) % Entre los voxs: 0.5404
[rho, pval] = corr(natfmu_aal, natfmu_aalMU) % Entre rois: 0.5083
[rho, pval] = corr(roinf, roinfMU) % Entre los voxs pero con valor de roi: 0.65 

figure, plot(natfmu)
hold on
plot(natfmu_Almu)

err = [natfci_aalMU(:,2) - natfci_aalMU(:,1)];
figure, errorbar(natfmu_aalMU, err, 'Marker', 'o','LineStyle','none', 'Color','k', 'MarkerFaceColor','b')
hold on
plot(natfmu_aal, 'Marker','s', 'MarkerFaceColor','r', 'LineStyle','none')
xticks(1:length(aal_label_reduc)), xticklabels(aal_label_reduc)
xlim([0 41]), xtickangle(90)

%% Mapa donde más diferencias por ROI hay

diffs = abs(natfmu_aal - natfmu_aalMU);

% Heatmap
figure('Position',[0 0 350 1000]), heatmap([natfmu_aal natfmu_aalMU], 'YDisplayLabels',aal_label_reducMU)
colormap("jet_omega_mod")
figure('Position',[300 0 10 1000]), heatmap(diffs, 'CellLabelColor','none')
colormap(flipud(colormap("gray")))
ax = gca; ax.ColorbarVisible = "off";

% Roi aissign
diffnf = zeros(1925,1);
for roi = 1:Nroi
    inds = find(label_inside_aal_reduc==roi);
    voxs = voxel_inside_aal(inds);
    diffnf(voxs) = diffs(roi);
end
% diffnf=max(diffnf)-diffnf % invertir para el colormap
Omega_nii(diffnf, 'Diffs_ROIs')

colcfg = [];
% colcfg.colim   = [0.7 3.4];  
colcfg.colmap   = flipud(colormap('gray'));
Omega_sourcefig(diffnf, 'test',colcfg)














%% Revisar lo que hay a partir de aquí











% Load roi str
load(['Z:\OMEGA\OMEGA-NaturalFrequencies-main\mat_files\aal_voxel_label_10mm.mat'])
Nroi = length(aal_label_reduc);

f1 = 0;
f2 = 3.7;
d  = 0.2;

foi1 = exp(f1:d:f2);
foi2 = exp(f1:d/10:f2);

for i = 1:length(foi2)
    c0 = foi2(i);
    rectp(i,:) = rectangularPulse(c0-c0/4 ,c0+c0/4 ,foi2);      % to isolate individual peaks later
end

p = parpool(4);                  % parallel pool
p.IdleTimeout = inf;

natf_aal = zeros(Nroi, size(natf.mu,2));
for it = 1:size(natf.mu,2)
for roi = 1:Nroi
    inds = find(label_inside_aal_reduc==roi);
    voxs = voxel_inside_aal(inds);
    fnatvox = natf.mu(voxs,it);
  
    [counts,~] = histcounts(fnatvox,foi1);
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
    % figure, hold on, plot(centers, counts)
    % plot(centers,squeeze(gausf),'k')

    i = find(rsq == max(rsq));      
    gausf = A(i) * exp(-(centers-mu(i)).^2 /(2*sigma(i).^2));
        % identify the one with the highest goodness of fit / amplitude
    natf_aal(roi,it)     = mu(i);
end
disp(['Iteration:' num2str(it) '/' num2str(size(natf.mu,2))])
end

roici(:,1) = prctile(natf_aal, 2.5,2)
roici(:,2) = prctile(natf_aal, 97.5,2)
roimu = median(natf_aal,2,"omitmissing");

roinf = zeros(1925,1);
for roi = 1:Nroi
    inds = find(label_inside_aal_reduc==roi);
    voxs = voxel_inside_aal(inds);
    roinf(voxs) = roimu(roi);
end

save roinatfreqgauss3 roimu roici natf_aal roinf


Omega_nii(roinf, 'Roi_Natfreq')

% CIs
figure, plot(roici)
difci = roici(:,2) - roici(:,1);
roici = zeros(1925,1);
for roi = 1:Nroi
    inds = find(label_inside_aal_reduc==roi);
    voxs = voxel_inside_aal(inds);
    roici(voxs) = difci(roi);
end
Omega_nii(roici, 'Roi_CI')



% Distribuciones
figure,
for rr = 1:Nroi
    subplot(4,10,rr)
% boxplot(natf_aal(rr,:))
histogram(natf_aal(rr,:))
title(aal_label_reduc{rr},'FontSize',8),
xlim([0 30])
end





















%% Natural frequencies (old)
freq_osc_time_cycles2 = freq_osc_time_secs;
freq_osc_time_cycles2(freq_osc_time_cycles2 == 0) = NaN;
all_natfreq_norm = normalize(freq_osc_time_cycles2,2,"zscore","robust");
avg_natfreq_norm = squeeze(mean(all_natfreq_norm,1,"omitmissing"));
nf= zeros(size(avg_natfreq_norm,1),1);

for v = 1:Nvoxin
    [pks,locs] = findpeaks(avg_natfreq_norm(v,:));
    if isempty(pks) | isinf(pks)
        nf(v) = NaN;
    else
        pks2 = pks(~isinf(pks));
        locs = locs(~isinf(pks));
        temp=frex(locs(find(pks2==max(pks2))));
        nf(v)=temp(1);
    end
end
% save natfeq.mat nf
Omega_nii(nf, 'Nf3cycmean')

% % Dominant frequencies
% avg_domfreq_norm = squeeze(mean(freq_osc_time_cycles,1));
% df = zeros(size(avg_domfreq_norm,1),1);
% 
% for v = 1:Nvoxin
%     [pks,locs] = findpeaks(avg_domfreq(v,:));
%     pks2 = pks(~isinf(pks));
%     locs = locs(~isinf(pks));
%     temp=frex(locs(find(pks2==max(pks2))));
%     df(v)=temp(1);
% end
% Omega_nii(nf, 'new_domfreqepis')

% Figure
colcfg.colim   = [0.7 3.4];  
colcfg.colmap   = 'jet_omega_mod';
Omega_sourcefig(log(natfmu), 'test',colcfg)


%% OLD
for roi = 1:Nroi
    inds = find(label_inside_aal_reduc==roi);
    voxs = voxel_inside_aal(inds);
    muvoxs = natf.mu(voxs,:);
    muvoxs(isnan(muvoxs)) =0; %%%%%%%%%% quito nans porque kmeans no converge
    
    Nk = 4;
    [idx,C,sumd,D] = kmeans(muvoxs,Nk,'Replicates',5,'MaxIter',200);      % k-means 4 clusters max
    nvoxs = [sum(idx==1) sum(idx==2) sum(idx==3) sum(idx==4) ];           % number of voxels in each cluster                       
    Nk = Nk - sum(nvoxs==1);                                              % eliminate clusters with only 1 member and repeat clustering                                           
    [idx,C,sumd,D] = kmeans(muvoxs,Nk,'Replicates',5,'MaxIter',200);      % k-means adapting number of clusters 
    nvoxs = [sum(idx==1) sum(idx==2) sum(idx==3) sum(idx==4) ];           % number of voxels in each cluster                       

    dominclust = find(nvoxs == max(nvoxs));                               % dominant cluster
    centroidvox = find(D(:,dominclust)==min(D(:,dominclust)));            % voxel closer to centroid
    represvox(roi) = voxs(centroidvox(1));                                % representative voxel for this AAL region
     % histogram(C'),title(aal_label_reduc{roi}),pause                    % visualize pattern of natural frequencies for both clusters
end

natfmu_aal      = median(freqvoxs(represvox,:),2, "omitmissing");
natfci_aal(:,1) = prctile(freqvoxs(represvox,:),2.5,2);       % 95% CI interval
natfci_aal(:,2) = prctile(freqvoxs(represvox,:),97.5,2);
% [natfmu_aal natfci_aal]

roinf = zeros(1925,1);
for roi = 1:Nroi
    inds = find(label_inside_aal_reduc==roi);
    voxs = voxel_inside_aal(inds);
    roinf(voxs) = natfmu_aal(roi);
end

Omega_nii(roinf, 'Roi_Natfreq')


% Correlación
idx = find(roinf>0);

corr(natfmu(idx), roinf(idx))
figure,plot(natfmu(idx),'.')
hold on
plot(roinf(idx),'.')

