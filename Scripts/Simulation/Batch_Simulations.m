%% Batch script to run simulations. 

%% 0. Prepare paths 
p.scripts = fileparts(fileparts(matlab.desktop.editor.getActiveFilename)); % Store the path of this script
addpath(genpath(p.scripts))

p.pBOSC = fileparts(p.scripts);
p.data = fullfile(p.pBOSC, 'Simresults');
p.files = fullfile(p.pBOSC, 'Files');

% Toolbox paths
addpath('Z:\Toolbox\fieldtrip-20230118') 
ft_defaults
addpath('Z:\Toolbox\')

%% 1. Input parameters

% Simulated voxels with oscillations
Nsources = Sim_voxels( [1700, 1567, 1409]);

% Simulated frequencies
Nfreqs =  [5 10 20];
% Simulated cycles
Ncycles = [3 10 20];
% Simulated SNR
Nsnr = [0.5 1 1.5];

repetitions = 5;

% Voxels with aperiodic component. 
% Semi-randomly select voxels at a distance from osicllatory voxels
% We generate 5 random voxels for each repetition (5 repetitions).
cfg = [];
cfg.numgens = 5;
cfg.repetitions = repetitions;
cfg.Nsources = Nsources;
cfg.seed = 211025; % seed for replication
evalc('apgens = Sim_apgenerators(cfg)');

%% Simulation loop
for vx = 1:length(Nsources) % Sources
    for fx = 1:length(Nfreqs) % Frequencies
        for cx = 1:length(Ncycles) % Cycles
            for snr = 1:length(Nsnr) % SNR
                for rep = 1:repetitions % repetitions
                    tic
                    fprintf([' \n Voxel: ' num2str(vx) '\n Frequency: ' num2str(fx) '\n Cycles: ' num2str(cx)  '\n SNR: ' num2str(snr) '\n Repetition: ' num2str(rep) '\n'])
                    cd(p.data)
                    mkdir(['Repetition_' num2str(rep)])
                    cd(['Repetition_' num2str(rep)])
                    %% Step 1: Simulate the signal
                    cfg = [];
                    cfg.dip = [Nsources(vx)];
                    cfg.apgens = apgens(rep,:);
                    cfg.simfreq = [Nfreqs(fx)];
                    cfg.simcycles = [Ncycles(cx)];
                    cfg.simsnr = [Nsnr(snr)];
               evalc('simsignal = pBOSC_SimulateSignal(cfg)');
                    %% Step 2: Source-reconstruction
               evalc('simsignal_source = pBOSC_SimulateBeamformer(simsignal)');
                    %% Step 3
                    cfg = [];
                    cfg.dip = [Nsources(vx)];
               evalc('sim_aperiodic = pBOSC_SimulateAperiodic(cfg, simsignal_source)');
                    %% Step 4
                    cfg = [];
                    cfg.dip = Nsources(vx);
                    cfg.ff = Nfreqs(fx);
                    cfg.cyc = Ncycles(cx);
                    cfg.snr = Nsnr(snr);
               evalc('[model] = pBOSC_SimulatedEpisodes(cfg, simsignal_source, sim_aperiodic)');

               % Save model performance
               save(['modelorig2_' num2str(vx), num2str(fx), num2str(cx), num2str(snr)], 'model');
                reptime = toc;
                fprintf([ ' \n Time of repetition: ' num2str( reptime / 60 ) ' min \n '])
                end
            end
        end
    end
end

