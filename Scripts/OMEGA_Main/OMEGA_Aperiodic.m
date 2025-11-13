function OMEGA_Aperiodic(sub, ses, dpath, winlen)

% sub: subect code, e.g., '0001'
% ses: session, e.g., '0001'
% dpath: folder to store processed data
% winlen = duration in seconds of the window to extract the aperiodic [Default is 20 seconds]

if ~exist('winlen')
    nsecs = 20;
else
    nsecs = winlen;
end

% -------------------------------------------------------- %
% 4.1 Fooof
% -------------------------------------------------------- %

% 0. Load data
load([dpath '\sub-' sub '\ses-' ses '\datasource_3423.mat'])

%% 4.1. FOOOF on 20 seconds length signal chunks

timelen = (dataclean.fsample * nsecs); % samples in nsecs (m)
windws = floor(length(dataclean.time{1}) / timelen);
t1 = 1;
aperiodic = [];

for tt = 1:windws
    disp(['Nb ' num2str(tt) ' of ' num2str(windws)])
    cfg = [];
    cfg.begsample = t1;

    t2 = t1 + timelen -1;
    cfg.endsample = t2;
    datachunk = ft_redefinetrial(cfg,dataclean);

    tmpdata = datachunk.trial{1};
    nfft = size(tmpdata,2);

    %% FFT
    data_fft = fft(tmpdata, [] ,2);

    if mod(nfft, 2) == 0
        amp_fft = abs(data_fft(:,1:nfft/2+1)) ./ nfft;
        frq = linspace(0, dataclean.fsample/2, nfft/2 + 1);
    else
        amp_fft = abs(data_fft(:,1:(nfft+1)/2)) ./ nfft;
        frq = linspace(0, dataclean.fsample/2, (nfft + 1)/2);
    end
    pow_fft = amp_fft.^2;

    if mod(nfft, 2) == 0
        amp_fft(:,2:end-1) = 2 * amp_fft(:,2:end-1);
        pow_fft(:,2:end-1) = 2 * pow_fft(:,2:end-1);
    else
        amp_fft(:,2:end) = 2 * amp_fft(:,2:end);
        pow_fft(:,2:end) = 2 * pow_fft(:,2:end);
    end

    %% 2. Fooof

    opts = getcfg_fooof(frq) % Get default brainstorm fooof cfg

    % Transform data format to apply fooof
    TF = zeros(size(pow_fft,1),1,size(pow_fft,2)-1); % Do not work with DC 0. Input all freqs but first
    pow_fft = pow_fft(:,2:end);
    TF(:,1,:) = double(movmean(pow_fft,10,2)); % smooth the powspctrm

    % fooof (from fieldtrip)
    [fs, fg] = process_fooof("FOOOF_matlab", TF, frq(2:end), opts, opts.hasOptimTools);

    % reconfigure structure
    apfit = [fg.ap_fit];
    apfit = reshape(apfit, [length(fg(1).ap_fit), length(fg)])';

    fooofedsp = [fg.fooofed_spectrum]; % model of powspctrm
    fooofedsp =  reshape(fooofedsp, [length(fg(1).fooofed_spectrum), length(fg)])';

    % peakfit = [fg.peak_fit]; % obtain it through foofedsp-apfit
    % peakfit = reshape(peakfit, [length(fg(1).peak_fit), length(fg)])';
    peakfit = fooofedsp - apfit;

    powspctmlog = [fg.power_spectrum]; % powspctrm en log10
    powspctmlog = reshape(powspctmlog, [length(fg(1).power_spectrum), length(fg)])';
    powspctm = 10.^powspctmlog; % original scale

    %% This powspctm matches exactly pow_fft, except the 1º value (DC)
    recap_pow = apfit;
    recap_pow = [pow_fft(:,1) recap_pow];
    recap_ampl= recap_pow;

    % Reconstruct 1/f^ time signal from the aperiodic fit
    if mod(nfft, 2) == 0
        recap_ampl(:,2:end-1) = recap_pow(:,2:end-1) ./ 2;
    else
        recap_ampl(:,2:end) = recap_pow(:,2:end) ./ 2;
    end

    recap_ampl = (sqrt(recap_ampl)) * nfft; %original ampl

    % Take the original signal phase
    phase = angle(data_fft(:,1:size(recap_ampl,2)));

    positivefx = recap_ampl .* exp(1i * phase);

    % Mirror positive frequencies (exlc. nyq)
    if mod(nfft, 2) == 0
        fullspectra = [ positivefx, conj(fliplr(positivefx(:,2:end-1)))];
    else
        fullspectra = [ positivefx, conj(fliplr(positivefx(:,2:end)))];
    end

    apreconsignal = real(ifft(fullspectra, nfft,2));
    aperiodic = [aperiodic apreconsignal];

end

save([dpath '\sub-' sub '\ses-' ses '\aperiodic.mat'], 'aperiodic')

end
%%

