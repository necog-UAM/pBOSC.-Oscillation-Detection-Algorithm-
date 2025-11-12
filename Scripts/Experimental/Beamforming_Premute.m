% cd('Z:\OMEGA\Enrique\Results')
% load dsimul.mat dataclean
cd('Z:\OMEGA\OMEGA_data\sub-0002\ses-0001')
load source_forward_10mm_allvox.mat
load dataclean.mat

grad  = dataclean.grad;
grid = source_forward.grid;
vol = source_forward.vol;

grad = ft_convert_units(grad, vol.unit);
grid = ft_convert_units(grid, vol.unit);

cfg             = [];
cfg.grid        = grid;
cfg.grad        = grad;
cfg.vol         = vol;
cfg.channel     = {'MEG'};
cfg.normalize   = 'no'; % <----
cfg.reducerank  = 2;
grid2           = ft_prepare_leadfield(cfg);

source_forward.grid = grid2;
% cd('Z:\OMEGA\Enrique\Results')
% save source_forward_sim source_forward

nperm = 500;
tic
for pp = 1:nperm
    permord = randperm(size(dataclean.trial{1},1));
    dataperm = dataclean;
    dataperm.trial{1} = dataperm.trial{1}(permord,:);
    

cfg            = [];
cfg.covariance = 'yes';
datacov        = ft_timelockanalysis(cfg, dataperm);       % covariance matrix

% Compute spatial filters (in source.avg.filter)
cfg                   = [];
cfg.method            = 'lcmv';
cfg.grad              = source_forward.grad;
cfg.headmodel         = source_forward.vol;
cfg.grid              = source_forward.grid;
cfg.lcmv.fixedori     = 'yes';
cfg.lcmv.normalize    = 'no'; % <----
cfg.lcmv.projectnoise = 'yes'; 
cfg.lcmv.keepfilter   = 'yes';          % important: save filters to use them later
cfg.lcmv.lambda       = '10%';          % the higher the smoother
cfg.lcmv.reducerank   = 2;
source                = ft_sourceanalysis(cfg, datacov);

load standard_sourcemodel3d10mm
source.avg.ori = {};
source.avg.mom = {};
source.avg.noisecov = {};
source.pos     = sourcemodel.pos;            % standard grid positions
source.inside  = source_forward.grid.inside;

% cd('Z:\OMEGA\Enrique\Results')
% save source_inverse_sim source

%% 3.5 Reconstruction of source-level activity
%% PERMUTE TO CenterHeadBias

time         = dataclean.time{1};
voxel_inside = find(source.inside==1);
Nvox         = length(voxel_inside);
datasourceperm   = zeros(Nvox,length(time));
datasource   = zeros(Nvox,length(time));
for i = 1:Nvox
    datasource(i,:) = source.avg.filter{voxel_inside(i)} * dataperm.trial{1};
end
datasourceperm = datasourceperm + datasource;
disp([num2str(pp) ' of ' num2str(nperm) ' permutations'])
end

datasourcepermmean = datasourceperm ./ nperm;
tperm = toc;
datasource = single(datasourcepermmean);

% Use only 5-minute recordings for all participants
if time(end) > 290      % >5 min
    t1 = findbin(time,290);
    t2 = findbin(time,time(end));
    time(t1:t2) = [];
    datasource(:,t1:t2) = [];
end

% Organize source-reconstructed data in a new fieldtrip structure
dataclean.trial{1}   = single(datasource);
dataclean.time{1}    = time;
dataclean.label      = {};
dataclean.sampleinfo = [1 size(dataclean.trial{1},2)];
for i = 1:Nvox
    dataclean.label{i} = ['V' num2str(i)];
end
clear datasource

cfg = [];
cfg.resamplefs = 256;
[dataclean] = ft_resampledata(cfg, dataclean)

cd('Z:\OMEGA\OMEGA_data\sub-0002\ses-0001')
    save 'datasource_perm' 'dataclean'
