%% Aperiodic test
% Sim signal
srate = 1000;
t = 1:1/srate:2.4;
% t = t(1:end-1);

sinewave = sin(2*pi*15*t) ;
figure, plot(t, sinewave)


pinknoise = f_alpha_gaussian(length(t), 1,1)';

signal = sinewave*0 + pinknoise;
figure, plot(signal)

%% Apply fooof

% FFT
    nfft = size(signal,2);

    data_fft = fft(signal, nfft,2) ./ nfft;

    frq = linspace(0, srate/2, floor(nfft/2)+1); %  nmber of freqs is nfft / 2 + 0 frq
    ampl_fft = abs(data_fft(:,1:length(frq))); % take the abs to get amplitude
    ampl_fft(:,2:end) = 2*ampl_fft(:,2:end); % double the positive excluding 0 dc
    powspctm_fft = ampl_fft.^2; % squared to get power

figure, plot(frq, powspctm_fft)
    TF = zeros(size(powspctm_fft,1),1,size(powspctm_fft,2)-1);

    ft_hastoolbox('brainstorm', 1);

    TF(:,1,:) = powspctm_fft(:,2:end);

    Freqs     = frq;
    Freqs(Freqs==0) = [];
    % This grabs the defaults from the brainstorm code
    opts_default  = getfield(process_fooof('GetDescription'), 'options');

    % Fetch user settings, this is a chunk of code copied over from
    % process_fooof, to bypass the whole database etc handling.
    cfg = [];
    opt                     = ft_getopt(cfg, 'fooof', []);
    opt.freq_range          = ft_getopt(opt, 'freq_range', frq([1 end]));
    opt.peak_width_limits   = ft_getopt(opt, 'peak_width_limits', opts_default.peakwidth.Value{1});
    opt.max_peaks           = ft_getopt(opt, 'max_peaks',         opts_default.maxpeaks.Value{1});
    opt.min_peak_height     = ft_getopt(opt, 'min_peak_height',   opts_default.minpeakheight.Value{1}/10); % convert from dB to B
    opt.aperiodic_mode      = ft_getopt(opt, 'aperiodic_mode',    opts_default.apermode.Value);
    opt.peak_threshold      = ft_getopt(opt, 'peak_threshold',    2);   % 2 std dev: parameter for interface simplification
    opt.return_spectrum     = ft_getopt(opt, 'return_spectrum',   1);   % SPM/FT: set to 1
    opt.border_threshold    = ft_getopt(opt, 'border_threshold',  1);   % 1 std dev: proximity to edge of spectrum, static in Python
    % Matlab-only options
    opt.power_line          = ft_getopt(opt, 'power_line',        'inf'); % for some reason it should be a string, if you don't want a notch, use 'inf'. Brainstorm's default is '60'
    opt.peak_type           = ft_getopt(opt, 'peak_type',         opts_default.peaktype.Value);
    opt.proximity_threshold = ft_getopt(opt, 'proximity_threshold', opts_default.proxthresh.Value{1});
    opt.guess_weight        = ft_getopt(opt, 'guess_weight',      opts_default.guessweight.Value);
    opt.thresh_after        = ft_getopt(opt, 'thresh_after',      true);   % Threshold after fitting always selected for Matlab (mirrors the Python FOOOF closest by removing peaks that do not satisfy a user's predetermined conditions)

    % Output options
    opt.sort_type  = opts_default.sorttype.Value;
    opt.sort_param = opts_default.sortparam.Value;
    opt.sort_bands = opts_default.sortbands.Value;

    % Check input frequency bounds
    if (any(opt.freq_range < 0) || opt.freq_range(1) >= opt.freq_range(2))
        bst_report('error','Invalid Frequency range');
        return
    end

    hasOptimTools = 0;
    if exist('fmincon', 'file')
        hasOptimTools = 1;
        disp('Using constrained optimization, Guess Weight ignored.')
    end

    [fs, fg] = process_fooof('FOOOF_matlab', TF, Freqs, opt, hasOptimTools);

    % reconfigure structure
    apfit = [fg.ap_fit];
    % apfit = reshape(apfit, [length(fg(1).ap_fit), length(fg)])';

    fooofedsp = [fg.fooofed_spectrum]; % model of powspctrm
    % fooofedsp =  reshape(fooofedsp, [length(fg(1).fooofed_spectrum), length(fg)])';

    % peakfit = [fg.peak_fit]; % obtain it through foofedsp-apfit
    % peakfit = reshape(peakfit, [length(fg(1).peak_fit), length(fg)])';
    peakfit = log10(fooofedsp) - log10(apfit);

    powspctm = [fg.power_spectrum]; % powspctrm en log10
    powspctm = reshape(powspctm, [length(fg(1).power_spectrum), length(fg)])';
    powspctm2 = 10.^powspctm; % original scale
%%

figure, plot(fs, powspctm,'k')
hold on, plot(fs, log10(apfit), 'r')
plot(fs, log10(fooofedsp),'b--')
xlim([0 35])
% figure, plot(fs, fooofedsp-apfit) % only peaks
% figure, plot(fs, log10(fooofedsp)-log10(apfit))

    % Reconstruct 1/f^ time signal
    tmp = sqrt(apfit);
    tmp = tmp ./2;
    tmp = [ampl_fft(1) tmp  tmp(end:-1:1)]; % zero for DC at the beginning.
    apreconsignal = bsxfun(@times, exp(bsxfun(@times,angle(data_fft),1i)), tmp);
    apreconsignal = ifft(apreconsignal,[],2,'symmetric') * nfft;

    % tmp = powspctm_fft;
    % tmp(1) = [];
    % tmp = sqrt(tmp);
    % tmp = tmp ./2;
    % tmp = [ampl_fft(1) tmp  tmp(end:-1:1)]; % zero for DC at the beginning.
    % apreconsignal = bsxfun(@times, exp(bsxfun(@times,angle(data_fft),1i)), tmp);
    % apreconsignal = ifft(apreconsignal,[],2,'symmetric') * nfft;
    
    figure, plot(t, signal)
    hold on
    plot(t, apreconsignal, '--o')

    aperiodic = [aperiodic apreconsignal];

   