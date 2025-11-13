function OMEGA_pBOSC(sub, ses, dpath)

% sub: subect code, e.g., '0001'
% ses: session, e.g., '0001'
% dpath: folder to store processed data

% -------------------------------------------------------- %
% 0  Load data
% 1  Prepare, baseline correct, downsample
% 2  Frequency analysis
% 3  Find power threshold
% 4  Find local peaks
% 5.5  Select detected peaks that exceed the power threshold
% 5.6  Find volume peaks
% 5.7  Create oscillatory episodes
% 5.8  Connect neighboring episodes
% -------------------------------------------------------- %

%% Step 0. Load data

% I. Load aperiodidc data to compute threshold
load([dpath '\sub-' sub '\ses-' ses '\aperiodic.mat'])

% Load original signal
load([dpath '\sub-' sub '\ses-' ses '\datasource_3423.mat'])


%% Step 1. Prepare and correct signal

% I. Cut original signal to same length of aperiodic
cfg = [];
cfg.begsample = 1;
cfg.endsample = length(aperiodic);
dataclean = ft_redefinetrial(cfg,dataclean);

% II. Correction to center of the head bias
dataclean.trial{1} = dataclean.trial{1} ./ rms(aperiodic,2);
noise_corr = aperiodic ./ rms(aperiodic,2);

dataclean_noise = dataclean;
dataclean_noise.trial{1} = noise_corr;

% III. Downsample both signals
cfg = [];
cfg.resamplefs = 128;
dataclean = ft_resampledata(cfg, dataclean);
dataclean_ap = ft_resampledata(cfg, dataclean_noise);
fsample = dataclean.fsample;

%% Step 2.Frequency analysis

% Frequency parameters
frex = exp(0.6:0.1:3.7); % 1.8 Hz to 40 Hz
Nfrex = length(frex);

% I. Original signal
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

% II. Noise signal
powspctmnoise = single(zeros(Nvox, Nfrex));
powspctmnoise = repmat(powspctmnoise, [1 1 Ntp]);

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
    freqnoise           = ft_freqanalysis(cfg, dataclean_ap);

    powspctmnoise(:,f0,:)= single(squeeze(freqnoise.powspctrm));
end

clear freq freqnoise

%% Step 3. Find the power threshold

% I. This thshld uses the 95% pctile of sim_aperiodic
thshld = single(zeros(Nfrex,1));
for f0 = 1:Nfrex
    thshld(f0) = prctile(prctile(squeeze(powspctmnoise(:,f0,:)),95,1),95,2);
end
thshld = repmat(thshld, 1, size(powspctmnoise,1));
thshld = permute(thshld,[2 1]);

save([dpath '\sub-' sub '\ses-' ses '\thshld.mat'], 'thshld')
clear powspctmnoise

%% Step 4. Find local peaks in each vox-freq point across time

[source] = Omega_source(size(powspctm,1), 1.5);
maxconn = max(sum(source.connmat));
detpks = int8(zeros(Nvox,Nfrex,Ntp));
%
for v = 1:Nvox
    vv = find(source.connmat(v,:)==1);
    if length(vv)==maxconn % maxconn=19
        [~,~,locmx] = findlocalmax(squeeze(powspctm(v,:,:)),1,[]);
        detpks(v,:,:) = locmx;
    end
end

%% Step 5. Select only peaks with power values above the APERIODIC threshold (95th percentile)

detpks_thap = int8(zeros(Nvox,Nfrex,Ntp));

for t=1:size(detpks,3)
    detpks2=detpks(:,:,t);
    powspctm2=powspctm(:,:,t);
    detpks2(powspctm2<thshld) = 0;
    detpks_thap(:,:,t)=detpks2;
end

%% Step 6. Select only tf points with a local maxima in the 3D brain volume

dim = source.dim;
xx = source.xx;
yy = source.yy;
zz = source.zz;

detpks2v = int8(zeros(1,Nfrex,Ntp));
detpks2v = repmat(detpks2v,[Nvox,1]);

for f0=1:Nfrex
    detpks2 = int8(zeros(dim(1),dim(2),dim(3),1));
    detpks2 = repmat(detpks2,[1,1,1,Ntp]);
    disp(['Frequency ' num2str(f0) '/' num2str(Nfrex)])

    dinterp = squeeze(powspctm(:,f0,:))'*source.dtempl;
    dinterp2 = reshape(dinterp,[Ntp dim(1),dim(2),dim(3)]);

    for t0=1:Ntp
        if sum(detpks(:,f0,t0))>0
            [pks,locs,locmx] = findlocalmax(squeeze(dinterp2(t0,:,:,:)), 3, []);
            [x,y,z] = ind2sub(size(locmx),find(locmx==1));
            for i = 1:length(pks)
                detpks2(x(i),y(i),z(i),t0) = 1;
            end
        end
    end

    detpks4D = int8(zeros(dim(1),dim(2),dim(3),1));
    detpks4D = repmat(detpks4D,[1,1,1,Ntp]);

    for v=1:Nvox
        detpks4D(xx(v),yy(v),zz(v),:) = detpks_thap(v,f0,:);
    end

    detpks2 = detpks2.*detpks4D;
    for v=1:Nvox
        detpks2v(v,f0,:) = detpks2(xx(v),yy(v),zz(v),:,:);
    end
end

Ntp = size(detpks2v,3);

clear detpks detpks4D detpks2

%% Step 7. Episode descriptives

[source1925] = Omega_source(1925, 1.5);

Nvoxin = length(source1925.inside);
voxin = find(ismember(source.inside,source1925.inside)); % The 1925 voxels inside the 3924

conn = conndef(2,'maximal');            % salen demasiados elementos conectados en el cluster 1; probar poniendo un thr arriba
epis = {};
timep_freq = (fsample).*(1./frex') ;    % time points for 1 cycles of each foi

for vin=1:Nvoxin   % only the 1925 voxels inside the cortex to save time
    disp(['Voxel ' num2str(vin) '/' num2str(Nvoxin)])
    v = voxin(vin);
    vv = find(source.connmat(v,:)==1);
    detpkstf = squeeze(logical(mean(detpks2v(vv,:,:),1)));    % use a sphere of 1.5 cm radius (smooth across voxels)  
    CC = bwconncomp(detpkstf,conn);
    L = labelmatrix(CC);

    ct=1;
    for i=1:CC.NumObjects
        if length(CC.PixelIdxList{i})>timep_freq(end)*1       % 1--- only evaluate clusters with duration > 3 cycles of highest freq
            temp = L==i;
            epis_freq = sum(temp,2);
            epis_freq = sum(epis_freq.*frex')./sum(epis_freq);          % weighted mean to get cluster frequency

            epis_tps = logical(sum(temp,1));
            epis_dur = sum(epis_tps);
            epis_dur_sec = epis_dur.*1./fsample;
            epis_dur_cyc = epis_dur_sec./(1./epis_freq);

            [f,t]=ind2sub(size(L),CC.PixelIdxList{i});
            v2 = repmat(v,[length(f),1]);
            if epis_dur_cyc >= 1                   % <----minimum duration of 1 cycle 
                epis{vin}(ct).freq = epis_freq;
                epis{vin}(ct).dur_sec = epis_dur_sec;
                epis{vin}(ct).dur_cyc = epis_dur_cyc;
                epis{vin}(ct).timeps = find(epis_tps==1);
                try
                    pow = powspctm(v2,f,t);
                catch
                    pow=[];
                    for j=1:length(v2)
                        pow(j)=powspctm(v2(j),f(j),t(j));
                    end
                end
                epis{vin}(ct).power = mean(pow(:));
                ct=ct+1;
            end
        end
    end
end

while length(epis)<Nvoxin
    epis{end+1}=[];
end

episocc = int8(zeros(Nvoxin,size(detpks2v,2),size(detpks2v,3)));      
for vin=1:Nvoxin   
    episocc(vin,:,:) = detpks2v(voxin(vin),:,:);
end

save([dpath '\sub-' sub '\ses-' ses '\episocc_matrix_rand.mat'], 'episocc') % matrix of episode occurrences

%% Step 8. Connect epis

cfg = [];
cfg.fsample = dataclean.fsample;
conepis = pBOSC_connect_epis(cfg, epis);

% Original epis: reduce volume and save in folder epis
epis2 = cell(size(epis));  % Preallocate cell array for epis2
for ep = 1:length(epis)
    fnams = fieldnames(epis{ep});
    epis2{ep} = repmat(struct(), size(epis{ep}));  % Preallocate structure array
    for sz = 1:length(epis{ep})
        for nams = 1:length(fnams)
            if strcmp(fnams{nams}, 'timeps')
                epis2{ep}(sz).(fnams{nams}) = int32(epis{ep}(sz).(fnams{nams}));
            else
                epis2{ep}(sz).(fnams{nams}) = single(epis{ep}(sz).(fnams{nams}));
            end
        end
    end
end

mkdir([dpath '\sub-' sub '\ses-' ses '\epis'])
cd([dpath '\sub-' sub '\ses-' ses '\epis'])

for vin = 1:Nvoxin
    epis = epis2{vin};
    filename = sprintf('epis_v%d.mat', vin);
    save(filename, 'epis');
end

% Connected epis: reduce volume and save in folder conepis
epis2 = cell(size(conepis));  % Preallocate cell array for epis2
for ep = 1:length(conepis)
    fnams = fieldnames(conepis{ep});
    epis2{ep} = repmat(struct(), size(conepis{ep}));  % Preallocate structure array
    for sz = 1:length(conepis{ep})
        for nams = 1:length(fnams)
            if strcmp(fnams{nams}, 'timeps')
                epis2{ep}(sz).(fnams{nams}) = int32(conepis{ep}(sz).(fnams{nams}));
            else
                epis2{ep}(sz).(fnams{nams}) = single(conepis{ep}(sz).(fnams{nams}));
            end
        end
    end
end

mkdir([dpath '\sub-' sub '\ses-' ses '\conepis'])
cd([dpath '\sub-' sub '\ses-' ses '\conepis'])

for vin = 1:Nvoxin
    epis = epis2{vin};
    filename = sprintf('epis_v%d.mat', vin);
    save(filename, 'epis');
end

% Only episodes > 3 cycles
epis3c = {};
Nvox = length(epis2);
for vox = 1:Nvox
    ct = 1;
    for i = 1:length(epis2{vox})
        if epis2{vox}(i).dur_cyc>=3                     %  only oscillatory episodes with > 3 cycles
            epis3c{vox}(ct).freq = epis2{vox}(i).freq;
            epis3c{vox}(ct).dur_sec = epis2{vox}(i).dur_sec;
            epis3c{vox}(ct).dur_cyc = epis2{vox}(i).dur_cyc;
            epis3c{vox}(ct).timeps = epis2{vox}(i).timeps;
            epis3c{vox}(ct).power = epis2{vox}(i).power;
            ct=ct+1;
        end
    end
end

cd([dpath '\sub-' sub '\ses-' ses])
save epis3cyc epis3c


clear powspctm detpks2v

end

