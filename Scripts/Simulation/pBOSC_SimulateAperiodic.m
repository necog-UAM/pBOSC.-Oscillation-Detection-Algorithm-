%% Step 3. Apply fooof in sliding windows and save the time aperiodic component

function sim_aperiodic = pBOSC_SimulateAperiodic(cfg, simsignal_source)

%--------------------------------------------------------%
% cfg.winlen = duration in seconds of the window to extract the aperiodic [Default is 20 seconds]
%--------------------------------------------------------%

% Default duration of the window is 20 seconds
if ~isfield(cfg, 'winlen')
    nsecs = 20;
else
    nsecs = cfg.winlen;
end
if ~isfield(cfg,'figflag')
    figflag = 0;
else
    figflag = cfg.figflag;
end
dip = cfg.dip;

%% 3.1. FOOOF on signal chunks

timelen = (simsignal_source.fsample * nsecs); % samples in nsecs (m)
windws = floor(length(simsignal_source.time{1}) / timelen);
t1 = 1;
aperiodic = [];

for tt = 1:windws
    tic
    cfg = [];
    cfg.begsample = t1;  
    t2 = t1 + timelen -1;
    cfg.endsample = t2;
    datachunk = ft_redefinetrial(cfg,simsignal_source);

    tmpdata = datachunk.trial{1};
    nfft = size(tmpdata,2);

    %% FFT
    data_fft = fft(tmpdata, [] ,2);

    if mod(nfft, 2) == 0
        amp_fft = abs(data_fft(:,1:nfft/2+1)) ./ nfft;
        frq = linspace(0, simsignal_source.fsample/2, nfft/2 + 1);
    else
        amp_fft = abs(data_fft(:,1:(nfft+1)/2)) ./ nfft;
        frq = linspace(0, simsignal_source.fsample/2, (nfft + 1)/2);
    end
    pow_fft = amp_fft.^2;

    if mod(nfft, 2) == 0
        amp_fft(:,2:end-1) = 2 * amp_fft(:,2:end-1);
        pow_fft(:,2:end-1) = 2 * pow_fft(:,2:end-1);
    else
        amp_fft(:,2:end) = 2 * amp_fft(:,2:end);
        pow_fft(:,2:end) = 2 * pow_fft(:,2:end);
    end

    % frq = linspace(0, datachunk.fsample/2, floor(nfft/2)+1); %  nmber of freqs is nfft / 2 + 0 frq

    %% 2. Fooof
   addpath('Z:\OMEGA\Enrique\OMEGA-Repo\Scripts\Basefunctions')
    opts = getcfg_fooof(frq) % Get default brainstorm fooof cfg

    % Transform data format to apply fooof
    TF = zeros(size(pow_fft,1),1,size(pow_fft,2)-1); % Do not work with DC 0. Input all freqs but first
    pow_fft = pow_fft(:,2:end);
    TF(:,1,:) = double(movmean(pow_fft,10,2)); % smooth the powspctrm

    % fooof
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
    % recap_pow = powspctm - peakfit;
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

    % Take the positive part of the signal with complex 0 phase
    % positivefx = recap_ampl .* exp(1i * zeros(size(recap_ampl)));
    % phase = 2 * pi * rand(size(recap_ampl)); % random
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
      
sim_aperiodic = aperiodic;

    % Figure: Aperiodic Power Spectrum
    if figflag
    figure, plot(sigampl(dip,2:end),'m', 'LineWidth',3, 'LineStyle','-.')
    % xlim([0 100]), ylim(1.0e+18*[0   1.2]),
    axis off
   
    % Figure: Aperiodic Time Signal
    figure, plot(apreconsignal(dip,:),'m'),
    axis off
    end

    
end


end


