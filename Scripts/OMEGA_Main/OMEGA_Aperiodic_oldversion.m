function OMEGA4_RemoveOscillations_ES (sub, ses, dpath)

% sub: subect code, e.g., '0001'
% ses: session, e.g., '0001'
% dpath: folder to store processed data

% sub: subect code, e.g., '0001'
% ses: session, e.g., '0001'
% dpath: folder to store processed data
% winlen = duration in seconds of the window to extract the aperiodic [Default is 20 seconds]
% Default duration of the window is 20 seconds

if ~exist('winlen')
    nsecs = 20;
else
    nsecs = winlen;
end

% -------------------------------------------------------- %
% 4.1 Fooof
% -------------------------------------------------------- %


% 0. Load data
cd([dpath '\sub-' sub '\ses-' ses])
sprintf(['\n \n' 'Processing Sub: ' num2str(sub) '...\n \n'])
load datasource_3423.mat

    %% 4.1. FOOOF on 20 seconds length signal chunks

    timelen = dataclean.fsample * nsecs; % samples in nsecs
    windws = floor(length(dataclean.time{1}) / timelen);
    t1 = 1;
    aperiodic = [];

    for tt = 1:windws
        tic
        disp(['Nb ' num2str(tt) ' of ' num2str(windws)])
        cfg = [];
        cfg.begsample = t1;
        t2 = t1 + timelen -1;
        t1 = t2 +1;
        cfg.endsample = t2;
        datachunk = ft_redefinetrial(cfg,dataclean);

        tmpdata = datachunk.trial{1};
        nfft = size(tmpdata,2);

        %% FFT
        data_fft = fft(tmpdata, nfft,2) ./ nfft;
           
        if mod(nfft, 2) == 0
            frq = linspace(0, dataclean.fsample/2, nfft/2 + 1);
            amp_fft = abs(data_fft(:,1:nfft/2+1));
        else
            amp_fft = abs(data_fft(:,1:(nfft+1)/2));
            frq = linspace(0, dataclean.fsample/2, (nfft + 1)/2);
        end
        if mod(nfft, 2) == 0
            amp_fft(:,2:end-1) = 2 * amp_fft(:,2:end-1);
            pow_fft = amp_fft.^2;
            TF = zeros(size(powspctm_fft,1),1,size(powspctm_fft,2)-1);
        else
            amp_fft(:,2:end) = 2 * amp_fft(:,2:end);
            pow_fft = amp_fft.^2;
            TF = zeros(size(powspctm_fft,1),1,size(powspctm_fft,2));
        end


        %% 2. Fooof
        opts = getcfg_fooof(frq); % Get default brainstorm fooof cfg

        TF(:,1,:) = double(movmean(pow_fft(:,2:end),10,2)); % smooth the powspctrm


        % fooof
        [fs, fg] = process_fooof("FOOOF_matlab", TF, frq(2:end), opts, opts.hasOptimTools);

        % reconfigure structure
        apfit = [fg.ap_fit];
        apfit = reshape(apfit, [length(fg(1).ap_fit), length(fg)])';

        fooofedsp = [fg.fooofed_spectrum];
        fooofedsp =  reshape(fooofedsp, [length(fg(1).fooofed_spectrum), length(fg)])';

        powspctm = [fg.power_spectrum];
        powspctm = reshape(powspctm, [length(fg(1).fooofed_spectrum), length(fg)])';
        powspctm2 = 10.^powspctm;

        % Reconstruct 1/f^ time signal
        tmp = sqrt(apfit);
        tmp = tmp ./2;
        tmp = [tmp zeros(size(tmp,1),1) tmp(:,end:-1:1)];
        apreconsignal = bsxfun(@times, exp(bsxfun(@times,angle(data_fft),1i)), tmp);
        apreconsignal = ifft(apreconsignal,[],2,'symmetric') * nfft;
        
        aperiodic = [aperiodic apreconsignal];
        toc
    end

    save('aperiodic_data', 'aperiodic')
end

