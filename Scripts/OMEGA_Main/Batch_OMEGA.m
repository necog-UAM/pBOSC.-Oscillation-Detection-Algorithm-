%% Batch script to run pBOSC in OMEGA data. 

%% ----------- 0. Prepare paths ---------------- %

p.scripts = fileparts(fileparts(matlab.desktop.editor.getActiveFilename)); % Store the path of this script
addpath(genpath(p.scripts))

p.pBOSC = fileparts(p.scripts); % Scripts folder
p.data = fullfile(p.pBOSC, 'OMEGA'); % Preprocessed OMEGA data folder
p.files = fullfile(p.pBOSC, 'Files'); % Essential files folder
p.raw = fullfile(p.pBOSC, 'Raw'); % Raw OMEGA data folder

% Toolbox paths
addpath('Z:\Toolbox\fieldtrip-20230118') % Path to Fieldtrip
ft_defaults

% Participant and session list
[Nsub, subs, sess] = Omega_subs();

%% ----------- 1. Reading, preprocessing and artifact correction ---------------- %

for s = 1:Nsub
    OMEGA_Preprocessing(subs{s}, sess{s}, p.raw, p.data);   
end

%% ----------- 2. Source reconstruction (computation of beamforming weights) ---------------- %

% MEG-MRI coregistration: mark fiducials (lpa, rpa and nasion; then quit) 
% and check result (if coregistration is not accurate, repeat procedure)
for s = 1:Nsub
    OMEGA_Coregistration(subs{s}, sess{s}, p.raw, p.data);         
end

% Compute forward model and beamforming weights
for s = 1:Nsub
    sprintf(['\n Processing Sub:' subs{s} '..\n'])
    evalc('OMEGA_Beamforming(subs{s}, sess{s}, p.data)');
end

%% ----------- 3. Get aperiodic to get thresholds ---------------- %

for s = 1:Nsub 
    sprintf(['\n Processing Sub:' subs{s} '..\n'])
    evalc('OMEGA_Aperiodic(subs{s}, sess{s}, p.data)');       
end

%% ----------- 4. pBOSC episode detection ---------------- %

for s = 1:Nsub
    sprintf(['\n Processing Sub:' subs{s} '..\n'])
    evalc('OMEGA_pBOSC(subs{s}, sess{s}, p.data)');   
end

