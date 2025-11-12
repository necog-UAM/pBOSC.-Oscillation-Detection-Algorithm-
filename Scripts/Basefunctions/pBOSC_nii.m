function pBOSC_nii(datasource, filename)

%--------------------------------------------%
% Generates a .nii file and saves it in current folder.
%--------------------------------------------%

[source] = Omega_source(size(datasource,1),1.5);

source.inverse.avg.pow(source.inside) = datasource;

source.inverse.avg.mom =cell(length(source.inverse.avg.pow),1);
source.inverse.time=1;

cfg = [];
cfg.parameter  = 'avg.pow';
cfg.downsample = 2;
cfg.interpmethod  =  'linear';
source_interp  = ft_sourceinterpolate (cfg, source.inverse, source.forward.mri);

cfg = [];
cfg.filetype  = 'nifti';
cfg.parameter = 'pow';
cfg.filename  = [filename];
ft_sourcewrite(cfg, source_interp)

end