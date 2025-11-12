function pBOSC_sourcefig(datasource, filename, varargin)

%--------------------------------------------%
% Generates a brain source image and saves it in current folder. Fieldtrip
% is needed. Can input color lims as cfg.colim.
%--------------------------------------------%

if nargin > 2
    optionalArg = varargin{1};
else
    optionalArg = 'default';
end

% Correct column orientation
dims = size(datasource);
if dims(2)>dims(1)
    datasource = datasource';
end

[source] = Omega_source(size(datasource,1),1.5);

source.inverse.avg.pow(source.inside) = datasource;
source.inverse.avg.mom = source.inverse.avg.filter;

cfg=[];
cfg.parameter  = 'avg.pow';
cfg.downsample = 2;
cfg.interpmethod = 'nearest';
source_interp = ft_sourceinterpolate (cfg, source.inverse, source.forward.mri);

figure('WindowState','maximized','Color',[1 1 1]);
% figure

cfg               = [];
cfg.method        = 'surface';
cfg.funparameter  = 'pow';
cfg.maskparameter = cfg.funparameter;
if strcmp(optionalArg, 'default')
    cfg.funcolorlim   = 'auto';
elseif isfield(optionalArg,'colim')
    cfg.funcolorlim   = optionalArg.colim;
end
if strcmp(optionalArg, 'default')
    cfg.funcolormap   = 'auto';
elseif isfield(optionalArg,'colmap')
    cfg.funcolormap   = optionalArg.colmap;
end
cfg.projmethod    = 'nearest';
cfg.opacity       = 0.8;
cfg.figure        = 'gca';
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

print('-dtiff','-r300',filename);

end