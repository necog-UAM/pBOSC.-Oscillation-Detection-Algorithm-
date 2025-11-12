clear all

addpath('Z:\OMEGA\Enrique\OMEGA-Repo\Scripts\Base Scripts')
rawpath = 'Z:\OMEGA\OMEGA_raw\';
dpath = 'Z:\OMEGA\OMEGA_data'

%% Step 0. Load data

    % I. Load empty-room or aperiodidc noise

    cd('Z:\OMEGA\OMEGA_data\sub-0002\ses-0001')
    load datasource_perm.mat
    dataperm = dataclean.trial{1};
    load datasource_allvox.mat
    load aperiodic.mat


     %% Step 1. Prepare and correct signal 
    
        % III. Correction to center of the head bias
        cfg = [];
        cfg.begsample = 1;
        cfg.endsample = length(aperiodic);
        dataclean = ft_redefinetrial(cfg,dataclean)

        dataperm = dataperm(:,1:size(dataclean.trial{1},2));
        dataclean.trial{1} = dataclean.trial{1} ./ rms(dataperm,2);
        noise_corr = aperiodic ./ rms(dataperm,2);

        % % % IV.a test with less time
        % sampwind = dataclean.fsample * 120; % in seconds
        % t1 = randi([1 length(dataclean.time{1})-sampwind]);
        % t2 = t1+sampwind;
        % cfg = [];
        % cfg.begsample = t1;
        % cfg.endsample = t2;
        % dataclean = ft_redefinetrial(cfg,dataclean)
        % dataclean_noise = dataclean;
        % dataclean_noise.trial{1} = noise_corr(:,t1:t2);

        % % Plot
        % figure, plot(dataclean.trial{1}(vxplot,:))
    
        % % IV.b Analyze Whole signal
        dataclean_noise = dataclean;
        dataclean_noise.trial{1} = noise_corr;

        % V. Downsample both signals
        cfg = [];
        cfg.resamplefs = 128;
        [dataclean] = ft_resampledata(cfg, dataclean);
        [dataclean_ap] = ft_resampledata(cfg, dataclean_noise);
        fsample = dataclean.fsample;
    
     %% Step 5.2.Frequency analysis
    
     % Load Frequencies
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
    
        % % Plot powspctm
        % vxplot = 741;
        % figure, imagesc(dataclean.time{1}, frex,squeeze(powspctm(vxplot,:,:))), title('Powspctm'), axis xy
        % figure, imagesc(dataclean_ap.time{1}, frex,squeeze(powspctmap(vxplot,:,:))), title('Aperiodic Powspctm'), axis xy
    
            clear freq freqnoise
    
     %% Step 1.2. Find the power threshold
  
     % This thshld uses the 95% pctile of aperiodic
     thshld = single(zeros(Nfrex,1));
     for f0 = 1:Nfrex
         thshld(f0) = prctile(prctile(squeeze(powspctmnoise(:,f0,:)),95,1),95,2);
     end

     thshld = repmat(thshld, 1, size(powspctmnoise,1));
     thshld = permute(thshld,[2 1]);

     save thshld thshld
    

        %% Step 1.3.Find local peaks in each vox-freq point across time

        [voxel_inside, Nvox, voxin, Nvoxin] = Omega_VoxInside();
        [dim,xx,yy,zz,connmat,dtempl,maxconn] = Omega_neighbors(1.5);

        detpks = int8(zeros(Nvox,Nfrex,Ntp));
        % 
        for v = 1:Nvox
            vv = find(connmat(v,:)==1);
            if length(vv)==maxconn%maxconn/2        %maxconn=19 nb %We divide /2 to allow vox with 10 neigh. To avoid spurious local max in corner voxels
                [~,~,locmx] = findlocalmax(squeeze(powspctm(v,:,:)),1,[]);
                detpks(v,:,:) = locmx;
            end
        end
    
        % Plot
        % figure,
        % vxplot = 2578;
        % clf,rtp = randi([1 size(detpks,3)], 1);
        % plot(1:length(frex),squeeze(powspctm(vxplot,:,rtp)))
        % [pks,locs,~] = findlocalmax(squeeze(powspctm(vxplpowspctmapot,:,:)),1,[]);
        % hold on,
        % plot(locs{rtp}, pks{rtp},'r*'),xticks(1:length(frex)), xticklabels(frex)
        % plot(thshld, 'k--')
    
        %% Step 1.4.Select only peaks with power values above the APERIODIC threshold (95th percentile)
    
        detpks_th = int8(zeros(Nvox,Nfrex,Ntp));

        for t=1:size(detpks,3)
            detpks2=detpks(:,:,t);
            powspctm2=powspctm(:,:,t);
            detpks2(powspctm2<thshld) = 0; 
            detpks_th(:,:,t)=detpks2;
        end
       
        clear thshld
        % Plot
        % figure, imagesc(dataclean.time{1},frex,squeeze(detpksth(vxplot,:,:))), axis xy, title('Detpks 2 pow threshold')
    
        %% Step 1.5.Select only tf points with a local maxima in the 3D brain volume
    
        detpks2v = int8(zeros(1,Nfrex,Ntp));
        detpks2v = repmat(detpks2v,[Nvox,1]);
    
        for f0=1:Nfrex
            detpks2 = int8(zeros(dim(1),dim(2),dim(3),1));
            detpks2 = repmat(detpks2,[1,1,1,Ntp]);
            disp(['Frequency ' num2str(f0) '/' num2str(Nfrex)])
         
            dinterp = squeeze(powspctm(:,f0,:))'*dtempl;
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
                detpks4D(xx(v),yy(v),zz(v),:) = detpks_th(v,f0,:);
            end
    
            detpks2 = detpks2.*detpks4D;
            for v=1:Nvox
                detpks2v(v,f0,:) = detpks2(xx(v),yy(v),zz(v),:,:);
            end
        end
    
        % Plot
        % figure, imagesc(dataclean.time{1},frex,squeeze(detpks2v(1794,:,:))), axis xy, title('Detpks 3 3D volume')
    
        clear detpks detpks4D detpks2
    
        %% Episode descriptives
        conn = conndef(2,'maximal');            % salen demasiados elementos conectados en el cluster 1; probar poniendo un thr arriba
        epis = {};
        timep_freq = (fsample).*(1./frex') ;    % time points for 1 cycles of each foi
    
        for vin=1:Nvoxin   % only the 1925 voxels inside the cortex to save time
            v = voxin(vin);
            vv = find(connmat(v,:)==1);
            detpkstf = squeeze(logical(mean(detpks2v(vv,:,:),1)));
            CC = bwconncomp(detpkstf,conn);
            L = labelmatrix(CC);
    
            % figure,imagesc(dataclean.time{1},log(frex),L),axis xy
            % set(gca,'YTickLabel',{round(exp(1:0.5:4))});
    
            ct=1;
            for i=1:CC.NumObjects
                if length(CC.PixelIdxList{i})>timep_freq(end)*1       %<3--- only evaluate clusters with duration > 3 cycles of highest freq
                    temp = L==i;
                    epis_freq = sum(temp,2);
                    epis_freq = sum(epis_freq.*frex')./sum(epis_freq);          % weighted mean to get cluster frequency
    
                    epis_tps = logical(sum(temp,1));
                    epis_dur = sum(epis_tps);
                    epis_dur_sec = epis_dur.*1./fsample;
                    epis_dur_cyc = epis_dur_sec./(1./epis_freq);
    
                    [f,t]=ind2sub(size(L),CC.PixelIdxList{i});
                    v2 = repmat(v,[length(f),1]);
                    if epis_dur_cyc >= 1                   % 3<----minimum duration of 3 cycles (see above)
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
            % figure,imagesc(dataclean.time{1},log(frex),L),axis xy
            % set(gca,'YTickLabel',{round(exp(1:0.5:4))});
        end
    
        while length(epis)<Nvoxin
            epis{end+1}=[];
        end
    
        %  % % Plot
        % epismatx = zeros(Nvoxin,Nfrex,Ntp);
        % for ii = 1:length(epis)
        %     tmptim = {};
        %     if isfield(epis{ii},'freq')
        %         tmpfrx = [epis{ii}.freq];
        %         for tt = 1:length(tmpfrx)
        %             idxfx = dsearchn(frex',tmpfrx')';
        %             epismatx(ii,idxfx(tt),[epis{ii}(tt).timeps]) = 1;
        %         end
        %     else
        %     end
        % end
        % vxplot=2425;
        % figure, imagesc(dataclean.time{1},frexA+-62º1º .0
         % ,squeeze(epismatx(dsearchn(voxin,vxplot),:,:)))
        % title(['Detpks vox ' num2str(vxplot) ' episodes']), axis xy
        % % 
        cd([dpath '\sub-' sub '\ses-' ses])
        % save('epis19', 'epis', '-v7.3') % para badsubs
        save epis19 epis
        clear powspctm detpks2v
        subtime = toc;
        disp(['Time for sub ' sub ' was ' num2str(subtime/3600) ' hours'])

