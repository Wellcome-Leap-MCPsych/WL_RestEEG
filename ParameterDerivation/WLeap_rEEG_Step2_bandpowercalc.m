%%% WellcomeLEAP - Resting State EEG data pre-processing %%%

% Last updated: 2024/01/18
% Author: Ty Lees

% This script requires:
% EEGLAB  and Matlab's Signal Processing Toolbox
% The following EEGLAB extensions/plugins (beyond the default)
% ERPLAB

% Use this code with pre-processed resting state data to:
% 1. Segment into appropriate windows
% 2. Calculate spectral band power in specified frequency bands
% 3. Export band power data

%% Workspace and variable set-up
clear variables;
eeglab          

cd('C:/Users/tlees/Partners HealthCare Dropbox/Ty Lees/WL_SMART_Trial/EEG Data Analyses/Baseline/eeg/EO Merged Data/');     %Set your CD to where you keep the data; this will work with cloud directories if you have the desktop apps
rawData = dir('*_processed.set');

%Pre-allocate arrays/tables that are iteratively filled
FilesProcessed = strings(length(rawData),1);

%% Load file and pre-process data
for subjID = 1:length(rawData)

    % Step 0: Load .set file
    loadName = strcat(rawData(subjID).folder,'\',rawData(subjID).name);  
    EEG = pop_loadset(loadName);                                            % Loads & reads processed EEG file

    % % Step 1: Events & Segmentation
    % % Step 1.1: Insert event marker every 2 seconds for segmenting
    % EEG = eeg_regepochs(EEG, 'recurrence', 2, 'rmbase', NaN, 'eventtype', 'Segment',...
    %                    'extractepochs', 'off');
    % 
    % % Step 1.2: Create-event list and segment data around inserted events
    % EEG = pop_epoch(EEG, {'Segment'}, [0  2], 'epochinfo', 'yes');
    %
    % % Step 1.3: Scan for and remove artifactual segments
    % EEG = pop_eegthresh(EEG, 1, 1:EEG.nbchan, -100, 100, EEG.xmin, EEG.xmax, 0, 1); % Removes epochs with voltages > ±100uV
    % EEG = pop_rejkurt(EEG, 1, 1:EEG.nbchan, 3, 3, 0, 1);    % Removes epochs with channels that have a kurtosis > ±3 SD

    % Step 2: Calculate band power
    % Step 2.1: Use Welch's method wih a 50% window overlap to calculate power up to Nyquist Freq in 0.5 Hz windows; 
    [spectra, freqs] = spectopo(EEG.data, 0, EEG.srate, 'winsize', EEG.srate*2, 'nfft', EEG.srate*2, 'overlap', EEG.srate, ...
                      'plot', 'off');   % Set second argument to 0 = use full data length; 
        
    % Step 2.2: Set frequency band sizes
    % Classic Bands
    deltaIdx = find(freqs>=1.5 & freqs<=6);
    thetaIdx = find(freqs >= 6.5 & freqs <= 8);
    alpha1Idx = find(freqs >= 8.5 & freqs <= 10);
    alpha2Idx = find(freqs >= 10.5 & freqs <= 12);
    beta1Idx = find(freqs >= 12.5 & freqs <= 18);
    beta2Idx = find(freqs >= 18.5 & freqs <= 21);
    beta3Idx = find(freqs >= 21.5 & freqs <= 30);
    gammaIdx = find(freqs >= 36.5 & freqs <= 44);

    % Alternative Bands
    deltaaltIdx = find(freqs >= 1.5 & freqs <= 4);
    thetaalt2Idx = find(freqs >= 4.5 & freqs <= 8);
    alphaaltIdx = find(freqs >= 8.5 & freqs <= 12);

    % Step 2.3 Calculate frequency band power per channel and convert from decibels to uV^2
    for i = 1:EEG.nbchan
        % Classic Bands
        deltaPower(:,i) = mean(10.^(spectra(i,deltaIdx)/10));
        thetaPower(:,i) = mean(10.^(spectra(i,thetaIdx)/10));
        alpha1Power(:,i) = mean(10.^(spectra(i,alpha1Idx)/10));
        alpha2Power(:,i) = mean(10.^(spectra(i,alpha2Idx)/10));
        beta1Power(:,i) = mean(10.^(spectra(i,beta1Idx)/10));
        beta2Power(:,i) = mean(10.^(spectra(i,beta2Idx)/10));
        beta3Power(:,i) = mean(10.^(spectra(i,beta3Idx)/10));
        gammaPower(:,i) = mean(10.^(spectra(i,gammaIdx)/10));

        % Alternate Bands
        deltaaltPower(:,i) = mean(10.^(spectra(i,deltaaltIdx)/10));
        thetaaltPower(:,i) = mean(10.^(spectra(i,thetaIdx)/10));
        alphaaltPower(:,i) = mean(10.^(spectra(i,alphaaltIdx)/10));
    end 
    
    % Step 3. Write spectral power output to files
    if not(isfolder([strcat(rawData(subjID).name(1:end-14), '_Output')]))
        mkdir([strcat(rawData(subjID).name(1:end-14), '_Output')])
    else 
        warning('Output directory already exists')
    end 
        
    % Classic Bands
    writematrix(deltaPower, strcat(cd, '\', rawData(subjID).name(1:end-14), '_Output\', rawData(subjID).name(1:end-4), '_DeltaPower.csv'));
    writematrix(thetaPower, strcat(cd, '\', rawData(subjID).name(1:end-14), '_Output\', rawData(subjID).name(1:end-4), '_ThetaPower.csv'));
    writematrix(alpha1Power, strcat(cd, '\', rawData(subjID).name(1:end-14), '_Output\', rawData(subjID).name(1:end-4), '_Alpha1Power.csv'));
    writematrix(alpha2Power, strcat(cd, '\', rawData(subjID).name(1:end-14), '_Output\', rawData(subjID).name(1:end-4), '_Alpha2Power.csv'));
    writematrix(beta1Power, strcat(cd, '\', rawData(subjID).name(1:end-14), '_Output\', rawData(subjID).name(1:end-4), '_Beta1Power.csv'));
    writematrix(beta2Power, strcat(cd, '\', rawData(subjID).name(1:end-14), '_Output\', rawData(subjID).name(1:end-4), '_Beta2Power.csv'));
    writematrix(beta3Power, strcat(cd, '\', rawData(subjID).name(1:end-14), '_Output\', rawData(subjID).name(1:end-4), '_Beta3Power.csv'));
    writematrix(gammaPower, strcat(cd, '\', rawData(subjID).name(1:end-14), '_Output\', rawData(subjID).name(1:end-4), '_GammaPower.csv'));

    % Alternate Bands
    writematrix(deltaaltPower, strcat(cd, '\', rawData(subjID).name(1:end-14), '_Output\', rawData(subjID).name(1:end-4), '_DeltaAltPower.csv'));
    writematrix(thetaaltPower, strcat(cd, '\', rawData(subjID).name(1:end-14), '_Output\', rawData(subjID).name(1:end-4), '_ThetaAltPower.csv'));
    writematrix(alphaaltPower, strcat(cd, '\', rawData(subjID).name(1:end-14), '_Output\', rawData(subjID).name(1:end-4), '_AlphaAltPower.csv'));
           
FilesProcessed(subjID,:) = EEG.setname; % Log the file that was just processed
end 

writematrix(FilesProcessed, strcat(cd,'\FilesProcessed_', string(date),'.txt')) %Writes list of files that were processed to a output file
