function opt = getcfg_fooof(frequency_range)

ft_hastoolbox('brainstorm', 1);

Freqs     = frequency_range;
Freqs(Freqs==0) = []; % eliminate the 0 dc

% Defaults from the brainstorm
defaultopts  = getfield(process_fooof('GetDescription'), 'options');

% Configure default and user options
opt                     = [];
opt.freq_range          = ft_getopt(opt, 'freq_range', frequency_range([1 end]));
opt.peak_width_limits   = ft_getopt(opt, 'peak_width_limits', defaultopts.peakwidth.Value{1});
opt.max_peaks           = ft_getopt(opt, 'max_peaks',         defaultopts.maxpeaks.Value{1});
opt.min_peak_height     = ft_getopt(opt, 'min_peak_height',   defaultopts.minpeakheight.Value{1}/10); % convert from dB to B
opt.aperiodic_mode      = ft_getopt(opt, 'aperiodic_mode',    defaultopts.apermode.Value);
opt.peak_threshold      = ft_getopt(opt, 'peak_threshold',    2);   % 2 std dev: parameter for interface simplification
opt.return_spectrum     = ft_getopt(opt, 'return_spectrum',   1);   % SPM/FT: set to 1
opt.border_threshold    = ft_getopt(opt, 'border_threshold',  1);   % 1 std dev: proximity to edge of spectrum, static in Python
% Matlab-only options
opt.power_line          = ft_getopt(opt, 'power_line',        'inf'); % for some reason it should be a string, if you don't want a notch, use 'inf'. Brainstorm's default is '60'
opt.peak_type           = ft_getopt(opt, 'peak_type',         defaultopts.peaktype.Value);
opt.proximity_threshold = ft_getopt(opt, 'proximity_threshold', defaultopts.proxthresh.Value{1});
opt.guess_weight        = ft_getopt(opt, 'guess_weight',      defaultopts.guessweight.Value);
opt.thresh_after        = ft_getopt(opt, 'thresh_after',      true);   % Threshold after fitting always selected for Matlab (mirrors the Python FOOOF closest by removing peaks that do not satisfy a user's predetermined conditions)

% Output options
opt.sort_type  = defaultopts.sorttype.Value;
opt.sort_param = defaultopts.sortparam.Value;
opt.sort_bands = defaultopts.sortbands.Value;

% Check input frequency bounds
if (any(opt.freq_range < 0) || opt.freq_range(1) >= opt.freq_range(2))
    bst_report('error','Invalid Frequency range');
    return
end

opt.hasOptimTools = 0;
if exist('fmincon', 'file')
    opt.hasOptimTools = 1;
    disp('Using constrained optimization, Guess Weight ignored.')
end

end