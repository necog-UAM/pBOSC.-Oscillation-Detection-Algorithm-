addpath(genpath('Z:\OMEGA\Enrique\OMEGA-Repo\Scripts\'))
addpath Z:\Toolbox\fieldtrip-20230118
ft_defaults
rawpath = 'Z:\OMEGA\OMEGA_raw\'; 
dpath = 'Z:\OMEGA\OMEGA_data\'
ft_warning off
addpath('Z:\Toolbox\NECOG')
addpath(genpath('Z:\OMEGA\OMEGA-NaturalFrequencies-main\mat_files'))

[Nsub, subs, sess] = Omega_subs(); %
[frex, Nfrex, fb_lab,fb1,fb2] = Omega_frex();
[source_forward, source,voxel_inside, dim,xx,yy,zz,connmat,dtempl] = Omega_voxin(1925, 1.5) ;
srate = 128;
Nvoxin = length(voxel_inside);

cd('Z:\OMEGA\Enrique\Results')
load oscillatory3_results_singlepis.mat 
[source_forward, source,voxel_inside, dim,xx,yy,zz,connmat,dtempl] = Omega_voxin(1925, 1.5) ;
[source_forward2, source2,voxel_inside2, ~,~,~,~,connmat2,dtempl2] = Omega_voxin(3423, 1.5) ;
Nvoxin = length(voxel_inside);
% Correspondence with voxels
for vx = 1:Nvoxin
correspidx(vx) = find(voxel_inside(vx)==voxel_inside2);
end

%% Matriz de los sujetos que más oscilan en cada frex
for fx = 1:32
    for s = 1:128
        [perc(s,fx),idx(s,fx)] = max(freq_osc_perc(s,:,fx));
    end
end
figure, imagesc(perc)
idx(43,8)
perc(43,8)

s=45;
cd(['Z:\OMEGA\OMEGA_data\sub-' subs{s} '\ses-' sess{s}])
load source_inverse_10mm_3423.mat
load datasource_3423vox.mat
inside = find(source.inside);
figure, plot(source.pos(inside,2), source.pos(inside,3),'ok', 'MarkerSize',10, 'LineWidth',2)

% figure, plot(dataclean.trial{1}(vxx,:))

% Save
vx = 43
hipvox = dataclean.trial{1}(correspidx(vx),:);
cd 'Z:\OMEGA\Enrique\Analysis Group'
save hipvox hipvox
%%
for fx = 1:32
    sub2most(:,fx) = squeeze(freq_osc_perc(3,:,fx));
end
figure, imagesc(sub2most)

tmp = squeeze(freq_osc_perc(2,:,3));
Omega_nii(tmp','s2')
frex(3)
%



ff=20;
frex(ff)
voxs = squeeze(freq_osc_perc(2,:,ff))';
Omega_nii(voxs, 'tst')

cd 'Z:\OMEGA\Enrique\Analysis Group'
parvox = dataclean.trial{1}(2693,:);
save parvox parvox

s=43;
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

s=s+1