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
        if isnan(max(single_frexvox(:,v)))
            freqvoxs(v,s) = NaN;
        else
            fc = find(single_frexvox(:,v) == max(single_frexvox(:,v))); % peak freq of the cluster with max Z-value
            if length(fc) > 1 % if more than one value, assign mode over neighbouring voxels
                vnb = find(connmat(v,:));
                for nb=1:length(vnb) % how many neighbors
                    tmpfc = find(single_frexvox(:,vnb(nb))==max(single_frexvox(:,vnb(nb))));
                    fcneigh(nb) = tmpfc(1); %%
                end
                fc = mode(fcneigh);
            end
            % if freq_osc_time_cycles(fc,v,s)==0 % If 0 total oscillations
            %     freqvoxs(v,s) = NaN;
            % end
            freqvoxs(v,s) = frex(fc);
        end
    end
end


%% Gaussian fit bootstrap

% foi2 = exp(0.6:0.02:3.7);
f1 = 0;
f2 = 3.7;
d  = 0.2;

foi1 = frex;
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
        nbvxs = connmat(v,:);
        fv = freqvoxs(nbvxs,randsub);
        [counts,~] = histcounts(fv(:),foi1);
        centers = exp(log(foi1) + d/2);  
        % clf,plot(centers(1:end-1),counts,'LineWidth',2), hold on

        counts  = interp1(centers(1:end-1),counts,foi2,'pchip');
        centers = foi2;
        % plot(centers,counts,'--')

        [pks,locs] = findpeaks(counts);
        [~,ii] = sort(pks,'descend');
        if length(ii) >= 3
            c0s = centers(locs(ii(1:3))); % fit a maximum of 3 peaks to speed up
        else
            c0s = centers(locs);
        end

        % gaussian fit of all candidate peak frequencies
        [A,mu,sigma,rsq] = gausfitc0(c0s,counts,centers,rectp,foi2);

        % i = find(rsq == max(rsq));% identify the one with the highest goodness of fit
        i = find(A == max(A));% identify the one with the highest goodness of fit
        natf.A(v,b)       = A(i);
        natf.mu(v,b)      = mu(i);
        natf.sigma(v,b)   = sigma(i);
        natf.rsq(v,b)     = rsq(i);
        natf.randsub(b,:) = randsub;
        % gausf = A(i) * exp(-(centers-mu(i)).^2 /(2*sigma(i).^2));
    % plot(centers, gausf,'k', 'LineWidth',2)
    % pause
    end
end

cd('Z:\OMEGA\Enrique\Results')
load natfreqgauss3_nb.mat
%       

natfmu      = median(natf.mu,2,"omitmissing");
natfstd = std(natf.mu,[],2, "omitmissing");
natfci(:,1) = prctile(natf.mu,2.5,2);          % 95% CI interval
natfci(:,2) = prctile(natf.mu,97.5,2);

Omega_nii(natfmu, 'nf_test_gauss2')

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


% Roi aissign
roinf = zeros(1925,1);
for roi = 1:Nroi
    inds = find(label_inside_aal_reduc==roi);
    voxs = voxel_inside_aal(inds);
    roinf(voxs) = natfmu_aal(roi);
end
Omega_nii(roinf, 'Experimental_ROI')

% LOAD Almu
load Z:\OMEGA\Enrique\Results\natfreq_aalMU.mat
natfmu_Almu      = median(natfMU.mu,2,"omitmissing");
natfciAlmu(:,1) = prctile(natfMU.mu,2.5,2);          % 95% CI interval
natfciAlmu(:,2) = prctile(natfMU.mu,97.5,2);

[rho,  pval] = corr(natfmu_Almu, natfmu) % Entre los voxs: 0.5404

colcfg = [];
colcfg.colim   = [0.7 3.4];  
colcfg.colmap   = 'jet_omega_mod';
Omega_sourcefig(log(natfmu), 'test',colcfg)



errAL = [natfci_aalMU(:,2) - natfci_aalMU(:,1)];
figure, errorbar(natfmu_aalMU, errAL/2, 'Marker', 'o','LineStyle','none', 'Color','b', 'MarkerFaceColor','b')
hold on
err = [natfci_aal(:,2) - natfci_aal(:,1)];
errorbar(natfmu_aal, err/2, 'Marker', 'o','LineStyle','none', 'Color','r', 'MarkerFaceColor','r')
xlim([0 41]),xticks(1:length(natfmu_aal)), xticklabels(aal_label_reducMU)
% If one contains the other
contains = [];
tol = 1e-4;
for vv = 1:length(natfmu_aal)
    v1 = natfci_aalMU(vv,1):.001:natfci_aalMU(vv,2);
    v2 = natfci_aal(vv,1):.001:natfci_aal(vv,2);
    intrsct = v1(ismembertol(v1,v2,tol));
    if isempty(intrsct)
        contains(vv) = NaN;
    else
        contains(vv)=vv;
    end
end
hold on
xline(contains,'Color', 'y','LineWidth',9,'Alpha',.2)
view([90 90])

%% Mapa donde más diferencias por ROI hay

diffs = abs(natfmu_aal - natfmu_aalMU);

% Heatmap
figure('Position',[0 0 350 1000]), heatmap([natfmu_aal natfmu_aalMU], 'YDisplayLabels',aal_label_reducMU)
colormap("jet_omega_mod")
figure('Position',[300 0 10 1000]), heatmap(diffs, 'CellLabelColor','none')
colormap(flipud(colormap("gray")))
ax = gca; ax.ColorbarVisible = "off";


