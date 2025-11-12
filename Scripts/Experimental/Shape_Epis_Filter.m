addpath(genpath('Z:\OMEGA\Enrique\OMEGA-Repo\Scripts'))

[Nsub, subs, sess] = Omega_subs(); %
[frex, Nfrex, fb_lab,fb1,fb2] = Omega_frex();
[source_forward, source,voxel_inside, dim,xx,yy,zz,connmat,dtempl] = Omega_voxin(1925, 1.5) ;
[source_forward2, source2,voxel_inside2, ~,~,~,~,connmat2,dtempl2] = Omega_voxin(3423, 1.5) ;
Nvoxin = length(voxel_inside);

s=2;
cd(['Z:\OMEGA\OMEGA_data\sub-' subs{s} '\ses-' sess{s}])
load singlepis.mat
load datasource_3423vox.mat
load('Z:\OMEGA\Enrique\Analysis Group\motorvox.mat')

% data = motorvox;
data = motorvox;
time = dataclean.time{1};
srate = dataclean.fsample;

% Correspondence with voxels
for vx = 1:Nvoxin
correspidx(vx) = find(voxel_inside(vx)==voxel_inside2);
end
% Correspondence in time
newtime = dataclean.time{1}(1:2:end);

%% Transform  epis from cells into a matrix

episodes = int8(zeros(Nvoxin,Nfrex,Ntp));
for v=1:Nvoxin
    for ep=1:length(epis{v})
        if epis{v}(ep).dur_cyc>=3 %<------------ 3 or 1 cycle episodes !!
            fm = epis{v}(ep).freq;
            fbin = dsearchn(frex',fm);
            tpts = epis{v}(ep).timeps;
            episodes(v,fbin,tpts) = 1;
        end
    end
end

vx = 1586;
vx3423 = correspidx(vx);

figure, surf(1:size(episodes,3), frex,squeeze(episodes(vx,:,:)))
view(0,90)
shading flat
axis tight

% Save timepoints where epis is detected as frex index
episidx = zeros(length(frex),length(time));
for ep = 1:length(epis{vx})
    for ffx = 1:length(frex)
    mainfx = dsearchn(frex',epis{vx}(ep).freq);
    if mainfx==ffx 
        tp = epis{vx}(ep).timeps;  
        sridx = tp(1)*2:tp(end)*2;
        if length(sridx) >= (3/mainfx*srate) % at least 3 cycle
            episidx(ffx,sridx) = frex(ffx);
        end
    end
    end
end
cd('Z:\OMEGA\Enrique\Analysis Group')
save motorepisidx episidx

% Save vox time
hipvox = dataclean.trial{1}(correspidx(vx),:);
cd 'Z:\OMEGA\Enrique\Analysis Group'
save hipvox hipvox


%% Animacion 1 second bandpass
timev = 1:128;
screen_size = get(0, 'ScreenSize');
figure('Position', [1, 1, screen_size(3), screen_size(4)]); 

for it = 1:200

    timecor = timev(2):2:timev(end);
    timecor = timecor/2;

    subplot(311)
    plot(timev, data(timev),'k'), xlim([timev(1) timev(end)])

    tmpidx = zeros(length(frex), length(timecor));
    for fx = 1:size(episodes,2)
        tmp = find(episodes(vx,fx,timecor));
        if ~isempty(tmp)
            tmpidx(fx,tmp')= 1;
            sridx = tmp(1)*2:tmp(end)*2;
            hold on, plot(timev(sridx),data(timev(sridx)))
        end
        fx
    end

    subplot(312)
    surf(1:length(timecor), frex,tmpidx)
    view(0,90)
    shading flat
    axis tight

    %% Filtering

    frx = find(sum(tmpidx,2));
    low_cutoff = frex(frx(1));   % Lower cutoff frequency (e.g., 100 Hz)
    high_cutoff =  frex(frx(end));  % Upper cutoff frequency (e.g., 300 Hz)

% Normalize the cutoff frequencies by the Nyquist frequency (Fs/2)
Wn = [low_cutoff high_cutoff] / (srate / 2);

% Create a Butterworth bandpass filter
[b, a] = butter(4, Wn, 'bandpass');  % 4th-order filter
subplot(313)
filtdata = filtfilt(b, a, data(timev));
plot(timev,filtdata,'LineWidth',2), xlim([timev(1) timev(end)])
title([num2str(low_cutoff) ' to ' num2str(high_cutoff) ' Bandpass'])


    timev = timev + 256;
    pause
clf
end


%% Filter by episodes

srate = dataclean.fsample;
figure
for ep = 1:length(epis{vx})
    tp = epis{vx}(ep).timeps;
    sridx = tp(1)*2:tp(end)*2;
    mainfx = dsearchn(frex',epis{vx}(ep).freq);
    if length(sridx) < (6/mainfx*srate)
        disp(num2str(ep))
    else
    % Filter
    low_cutoff = frex(mainfx-1);  
    high_cutoff =  frex(mainfx)+1;  

    Wn = [low_cutoff high_cutoff] / (srate / 2);
    [b, a] = butter(4, Wn, 'bandpass');  % 4th-order filter

    % % zeropad
    % pad_len = 100;
    % x_padded = [zeros(1,pad_len) data(sridx) zeros(1,pad_len)];
    filtdata = filtfilt(b, a, data(sridx));
    % filtdata = filtdata(pad_len+1:end-pad_len); % Remove padding

    plot(data(sridx),'k')
    hold on
    plot(filtdata,'LineWidth',2)
    title(num2str([epis{vx}(ep).freq]))
    pause
    clf
    end
end

