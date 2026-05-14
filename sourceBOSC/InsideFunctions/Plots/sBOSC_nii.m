function sBOSC_nii(datasource, filename, cfg)
%%
% Generates a .nii file from source data and saves it in the current folder.
% Fieldtrip required.
%
%   sBOSC_nii(datasource, filename) generates and saves the NIfTI file.
%   sBOSC_nii(datasource, filename, cfg) applies custom configurations.
%
%   Inputs:
%           datasource - Numeric vector containing exactly 1925 or 3294 voxels.
%           filename   - String or char array specifying the output .nii file name.
%           cfg        - (Optional) Structure containing configurations:
%               cfg.interp     : String defining interpolation method. Default: 'nearest'.
%               cfg.downsample : Numeric scalar for downsampling. Default: 2.

% Fieldtrip path check
ftPath = fileparts(which('ft_defaults'));
if ~isempty(ftPath)
    addpath(genpath(ftPath));
end

% Defaults
if nargin < 2 || isempty(filename)
    error('A filename is required to save a NIfTI file.');
end
if nargin < 3
    cfg = struct();
end
if isfield(cfg, 'interp')
    interpmethod = cfg.interp;
else
    interpmethod = 'nearest';
end
if isfield(cfg, 'downsample')
    ds_factor = cfg.downsample;
else
    ds_factor = 2;
end

% Correct column orientation
datasource = datasource(:);
Nvox = length(datasource);

if Nvox ~= 1925 && Nvox ~= 3294
    error('Input data must have exactly 1925 or 3294 voxels.');
end

% Load sources
if Nvox == 1925
    t = load('source_template_10mm_1925.mat');
else
    t = load('source_template_10mm_3294.mat');
end
inside = t.source.inverse.inside;

mdl = load('standard_sourcemodel3d10mm.mat');
source_mni = mdl.sourcemodel;
source_mni.inside = inside;
source_mni.avg.pow = nan(length(inside), 1);
source_mni.avg.pow(inside) = datasource;
source_mni.time = 1;

templ = load('standard_mri.mat');

% Interpolation
cfg_interp = [];
cfg_interp.parameter    = 'avg.pow';
cfg_interp.downsample   = ds_factor;
cfg_interp.interpmethod = interpmethod;
source_interp = ft_sourceinterpolate(cfg_interp, source_mni, templ.mri);

% Write to NIfTI
cfg_write = [];
cfg_write.filetype  = 'nifti';
cfg_write.parameter = 'pow';
cfg_write.filename  = filename;
ft_sourcewrite(cfg_write, source_interp);
end