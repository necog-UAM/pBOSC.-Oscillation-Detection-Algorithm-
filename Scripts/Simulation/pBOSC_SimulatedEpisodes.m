%% Step 4. Frequency decomposition and episode detection

function [model] = pBOSC_SimulatedEpisodes(cfg, simsignal_source, sim_aperiodic)

%
if ~isfield(cfg,'figflag')
    figflag = 0;
else
    figflag = cfg.figflag;
end
%
 dip = cfg.dip;
 ff = cfg.ff;
 cyc = cfg.cyc;
 snr = cfg.snr;
%

% Frequency parameters
frex = exp(0.6:0.1:3.7); % 1.8 Hz to 40 Hz
Nfrex = length(frex);

source = Omega_source(size(simsignal_source.trial{1},1),1.5);
inside = find(source.inverse.inside);

%% Step 0. Prepare and correct signal
if size(simsignal_source.trial{1},2)>size(sim_aperiodic,2)
    % I. Cut original signal to same length of sim_aperiodic
    cfg = [];
    cfg.begsample = 1;
    cfg.endsample = size(sim_aperiodic,2);
    simsignal_source = ft_redefinetrial(cfg,simsignal_source)
else
    sim_aperiodic = sim_aperiodic(:,1:size(simsignal_source.trial{1},2));
end
    % II. Correction to center of the head bias
    simsignal_source.trial{1} = simsignal_source.trial{1} ./ rms(sim_aperiodic,2);
    noise_corr = sim_aperiodic ./ rms(sim_aperiodic,2);
    
    % III. Analyze Whole signal
    simsignal_source_noise = simsignal_source;
    simsignal_source_noise.trial{1} = noise_corr;
    
    % IV. Downsample both signals
    cfg = [];
    cfg.resamplefs = 128;
    [simsignal_source] = ft_resampledata(cfg, simsignal_source);
    [simsignal_source_ap] = ft_resampledata(cfg, simsignal_source_noise);
    fsample = simsignal_source.fsample;

%% Step 1.Frequency analysis
    
     % I. Original signal
        Nvox = length(simsignal_source.label);
        Ntp = length(simsignal_source.time{1});
        powspctm = single(zeros(Nvox, Nfrex));
        powspctm = repmat(powspctm, [1 1 Ntp]);
    
        for f0 = 1:Nfrex
            cfg            = [];
            cfg.method     = 'mtmconvol';
            cfg.taper      = 'hanning';
            cfg.output     = 'pow';
            cfg.foi        = frex(f0);
            cfg.toi        = simsignal_source.time{1};
            cfg.t_ftimwin  = 5./frex(f0);
            cfg.pad        = 'nextpow2';
            cfg.keeptrials = 'yes';
            freq           = ft_freqanalysis(cfg, simsignal_source);
    
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
            cfg.toi        = simsignal_source_ap.time{1};
            cfg.t_ftimwin  = 5./frex(f0);
            cfg.pad        = 'nextpow2';
            cfg.keeptrials = 'yes';
            freqnoise           = ft_freqanalysis(cfg, simsignal_source_ap);
    
            powspctmnoise(:,f0,:)= single(squeeze(freqnoise.powspctrm));
        end
        
            clear freq freqnoise

            % Figure: Power spectrum
            if figflag
            figure, imagesc(squeeze(powspctm(dip,:,:))), colormap('sky'), axis xy,
            cx = caxis(), axis off
            figure, imagesc(squeeze(powspctmnoise(dip,:,:))), colormap('sky'), axis xy,
            caxis(cx), axis off

            figure('Position',[0 0 1000 800])
            plot(source.inverse.pos(inside,2), source.inverse.pos(inside,3),'ok', 'MarkerSize',15, 'LineWidth',2)
            set(gca, 'Ylim', [-12 14])
            set(gca, 'Xlim', [-14 10])
            axis off
            end

%% Step 2. Recreate the original simulated signal
    simepisodes = int8(zeros(Nvox,Nfrex,Ntp));

    dipolevox = simsignal_source.sim.dipole;
    dipolefrex = simsignal_source.sim.freq;
    dipolecyc = simsignal_source.sim.cycles;
        dwnsmpl = simsignal_source.sim.fsample ./ simsignal_source.fsample;
    dipolepnts = downsample(simsignal_source.sim.osctimepoints, dwnsmpl);
    dipolepnts = floor(dipolepnts/dwnsmpl);

    simepisodes(dipolevox, dsearchn(frex', dipolefrex), dipolepnts) = 1;

%% Step 3. Find the power threshold
  
     % I. This thshld uses the 95% pctile of sim_aperiodic
     thshld = single(zeros(Nfrex,1));
     for f0 = 1:Nfrex
         thshld(f0) = prctile(prctile(squeeze(powspctmnoise(:,f0,:)),95,1),95,2);
     end
     thshld = repmat(thshld, 1, size(powspctmnoise,1));
     thshld = permute(thshld,[2 1]);

     % Figure: Aperiodic power threshold
     if figflag
     f0=18;
     c0=dip;
     datadist = squeeze(powspctmnoise(c0,f0,:));
     perc95 = prctile(datadist, 95);
     h =  histogram(datadist)  ;
     h.EdgeColor = 'm'
     h.FaceColor = 'm'
     h.FaceAlpha = .3
     hold on, xline(perc95,'k--', 'LineWidth',5)
     axis off
     figure, plot(simsignal_source.trial{1}(dip,:),'k', 'LineWidth',1), axis off
     end
 
%% Step 4. Find local peaks in each vox-freq point across time 

        detpks = int8(zeros(Nvox,Nfrex,Ntp));
        maxconn = max(sum(source.connmat));
        for v = 1:Nvox
            vv = find(source.connmat(v,:)==1);
            if length(vv)==maxconn % To avoid spurious local max in corner voxels
                [pkk,lkk,locmx] = findlocalmax(squeeze(powspctm(v,:,:)),1,[]);
                detpks(v,:,:) = locmx;
            end
        end
  
%% Step 5. Select only peaks with power values above the sim_aperiodic threshold (95th percentile)
    
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
        
        clear detpks detpks4D detpks2

%% Step 7. Episode descriptives (only 1925 inside voxels).
        conn = conndef(2,'maximal'); 
        epis = {};
        timep_freq = (fsample).*(1./frex');    % time points for 1 cycles of each foi

[source1925] = Omega_source(1925,1.5);

        Nvoxin = length(source1925.inside);

       voxin = find(ismember(source.inside,source1925.inside)); % The 1925 voxels inside the 3924
       for vin=1:Nvoxin   % only the 1925 voxels inside the cortex to save time
            v = voxin(vin);
            vv = find(source.connmat(v,:)==1);
            detpkstf = squeeze(detpks2v(v,:,:));
            CC = bwconncomp(detpkstf,conn);
            L = labelmatrix(CC);
        
            ct=1;
            for i=1:CC.NumObjects
                if length(CC.PixelIdxList{i})>timep_freq(end)*1 %< only evaluate clusters with duration > 1 cycle of highest freq
                    temp = L==i;
                    epis_freq = sum(temp,2);
                    [~,idx_freq] = max(epis_freq);
                    epis_freq = frex(idx_freq);
                    epis_tps = logical(sum(temp,1));
                    epis_dur = sum(epis_tps);
                    epis_dur_sec = epis_dur.*1./fsample;
                    epis_dur_cyc = epis_dur_sec./(1./epis_freq);
    
                    [f,t]=ind2sub(size(L),CC.PixelIdxList{i});
                    v2 = repmat(v,[length(f),1]);
                    if epis_dur_cyc >= 1                   % minimum duration of 1 cycle
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
                i=i+1;
            end
        end
    
        while length(epis)<Nvoxin
            epis{end+1}=[];
        end

            % Transform episodes from cells back to matrix

    episodes = int8(zeros(Nvoxin,Nfrex,Ntp));
    for v=1:Nvoxin
        for ep=1:length(epis{v})
            fm = epis{v}(ep).freq;
            fbin = dsearchn(frex',fm);
            tpts = epis{v}(ep).timeps;
            episodes(v,fbin,tpts) = 1;
        end
    end

    % Transform simulated episodes to inside voxels to match episodes size
    simepisodesIN = simepisodes(voxin,:,:);
    dipolevoxIN = find(voxin==dipolevox); % Voxel with osillation in the 1925 inside voxs

%% Step 8. Connect epis

cfg = [];
cfg.fsample = fsample;
conepis = pBOSC_connect_epis(cfg, epis);


    % Transform episodes from cells back to matrix

    conepisodes = int8(zeros(Nvoxin,Nfrex,Ntp));
    for v=1:Nvoxin
        for ep=1:length(conepis{v})
            fm = conepis{v}(ep).freq;
            fbin = dsearchn(frex',fm);
            tpts = conepis{v}(ep).timeps;
            conepisodes(v,fbin,tpts) = 1;
        end
    end

%% Step 9. Compare results
  conepismodel = pBOSC_Stats(conepisodes, simepisodesIN, powspctm(voxin,:,:), dipolevoxIN, figflag);

  % Check real distance between predicted and simulated
  simvx = str2double(conepismodel.simvoxel{8,2});
  altvx = str2double(conepismodel.altvoxel{8,2});

  pos = source1925.inverse.pos(source1925.inside,:);

  % Calculate distance between simulated voxel and detected voxel
  altDist = pdist2(pos(simvx,:), pos(altvx,:));

  conepismodel.altvoxel{9,1} = ['Alt voxel Distance: '];
  conepismodel.altvoxel{9,2} = [num2str(altDist)];
  
    model = conepismodel;

end