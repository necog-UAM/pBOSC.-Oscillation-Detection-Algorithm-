addpath(genpath('Z:\OMEGA\Enrique\OMEGA-Repo\Scripts\'))
rawpath = 'Z:\OMEGA\OMEGA_raw\'; 
dpath = 'Z:\OMEGA\OMEGA_data\'

[Nsub, subs, sess] = Omega_subs(); %
[frex, Nfrex, fb_lab,fb1,fb2] = Omega_frex();
[source_forward, source,voxel_inside, dim,xx,yy,zz,connmat,dtempl] = Omega_voxin(1925, 1.5) ;
srate = 128;
Nvoxin = length(voxel_inside);

% cd('Z:\EEG_DATA\OMEGA\Enrique\Results')
load oscillatory3_results_singlepis.mat 

%% Matriz de los sujetos que más oscilan en cada frex
for fx = 1:32
    for s = 1:128
        [perc(s,fx),idx(s,fx)] = max(freq_osc_perc(s,:,fx));
    end
end
figure, imagesc(perc)
idx(103,25)
perc(57,28)

s=2;
cd(['Z:\OMEGA\OMEGA_data\sub-' subs{s} '\ses-' sess{s}])
load source_inverse_10mm_3423.mat
load datasource_3423vox.mat
inside = find(source.inside);
figure, plot(source.pos(inside,2), source.pos(inside,3),'ok', 'MarkerSize',10, 'LineWidth',2)

% figure, plot(dataclean.trial{1}(vxx,:))

% Save
frontvox = dataclean.trial{1}(475,:);
cd 'Z:\OMEGA\Enrique\Analysis Group'
save hipvox hipvox







ff=20;
frex(ff)
voxs = squeeze(freq_osc_perc(2,:,ff))';
Omega_nii(voxs, 'tst')

%% Save episodes
for s=1:Nsub;
    cd(['Z:\OMEGA\OMEGA_data\sub-' subs{s} '\ses-' sess{s}])
    load singlepis.mat
    episodes = int8(zeros(Nvoxin,Nfrex,Ntp));
    for v=1:Nvoxin
        for ep=1:length(epis{v})
            if epis{v}(ep).dur_cyc>=1 %<------------ 3 or 1 cycle episodes !!
                fm = epis{v}(ep).freq;
                fbin = dsearchn(frex',fm);
                tpts = epis{v}(ep).timeps;
                episodes(v,fbin,tpts) = 1;
            end
        end
    end
    save epismatx episodes
    s=s+1
end
