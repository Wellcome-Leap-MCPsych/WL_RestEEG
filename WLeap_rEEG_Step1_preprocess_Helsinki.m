%%% WellcomeLEAP - Resting State EEG data pre-processing %%%
%%% For Helsinki

% Last updated: 2024/01/16
% Author: Ty Lees

% This script requires:
% EEGLAB 
% The following extensions/plugins (beyond the default)
% dipfit, firfilt, ERPLAB, FASTER, ICLabel, bva-io

% Use this code with raw recorded data to:
% 1. Import data into EEGLAB
% 2. Check for and/or add channel information
% 3. Optional resampling - currently checks for a SR of 1000
% 4. Bad Channel rejection using FASTER
% 5. Bandpass & Notch Filtering
% 6. Re-reference the data
% 7. ICA decomposition, labelling, and removal
% 8. Interpolate & zero the removed bad channels

%% Workspace and variable set-up
clear variables;
eeglab          

cd('...');     %Set your CD to where you keep the data; this will work with cloud directories if you have the desktop apps
rawData = dir('*.vhdr');

%Pre-allocate arrays/tables that are iteratively filled
FilesProcessed = strings(length(rawData),1);

%% Load file and pre-process data
for subjID = 1:length(rawData) 

    % Step 1: Import data; function depends on data format
    EEG = pop_loadbv(rawData(subjID).folder, rawData(subjID).name); % Load data for Helsinki, Stanford, McLean
    EEG.setname = rawData(subjID).name(1:end-5);    % Create .set name for display in EEGLab
    EEGPreprocessTAB = table({EEG.setname},'VariableNames',{'Participant'});            % Record participant ID and create table for processing data
    EEG = eeg_checkset(EEG); % Check dataset integrity

    % Step 2: Add Channel information and remove unwanted channels  
    EEG = pop_chanedit(EEG, 'load',{strcat(cd,'\WL_Helsinki_ChanLoc_Info.ced'),'filetype','autodetect'});  % Helsinki 
    EEG = pop_select(EEG,'nochannel',{'EOG1','EOG2'}');   

    % Step 3: Resample to 1000 Hz if needed
    EEGPreprocessTAB.Orig_srate = EEG.srate;    % Record original sample rate

    if EEG.srate > 1000          
        EEG = pop_resample(EEG, 1000);
    elseif EEG.srate < 1000
        warning('Sample rate below target; function is trying to up-sample')
        fid = fopen(strcat(rawData(subjID).name,'_Resample_warning.txt'), 'wt');
            fprintf(fid, 'Recording sample rate is below target, and function is trying to up-sample');
        fclose(fid);
    end
    
    % Step 4: Remove bad channels - FASTER method - Requires all channels to have location info
    EEG = pop_reref(EEG, {'Pz'}, 'keepref','on');       % Temporarily re-reference to Pz; required for FASTER

    list_properties = channel_properties(EEG, 1:EEG.nbchan, 49);     % Computes bad channel coefficents; Arguments are: Data input, Which Chans, Ref Chan ID   
    badChanindices = min_z(list_properties);
    badChans = find(badChanindices == 1);   % Identifies bad channels
    badChans = badChans(badChans ~= 49);     % Removes the reference from list so it's not excluded   					
    
    EEGPreprocessTAB.num_bad_chans = numel(badChans);       % Record number of bad channels identified
    EEGPreprocessTAB.badchan_indices = num2str(reshape(badChans',1,[]));    % Record index of each bad channel
    
    EEG = pop_select(EEG,'nochannel', badChans);    % Remove bad channels from recording    

    % Step 5: Filter data
        % Step 5.1: Bandpass filter at 0.1 - 80 Hz
        EEG  = pop_basicfilter(EEG, 1:EEG.nbchan, 'Boundary', 'boundary', 'Cutoff', [ 0.1 80], 'Design', 'butter', 'Filter', 'bandpass', 'Order', 2, 'RemoveDC', 'on');
        % Step 5.2: Notch filter @ 50 Hz
        EEG  = pop_basicfilter(EEG, 1:EEG.nbchan, 'Boundary', 'boundary', 'Cutoff',  50, 'Design', 'notch', 'Filter', 'PMnotch', 'Order', 180); 

    % Step 6: Re-reference to common average
    EEG = pop_reref(EEG, []);

    % Step 7: ICA
        % Step 7.1: Run decomposition using FastICA.m
        EEG = pop_runica(EEG,'icatype','fastica','approach','symm','g','tanh','stabilization','on');
    
        % Step 7.2: Use ICLabel to label components, flag artifactual componens, collate information
        EEG = pop_iclabel(EEG, 'Default'); 
        
        EEG = pop_icflag(EEG, [NaN NaN;0.9 1;0.9 1;NaN NaN;0.9 1;0.9 1;NaN NaN]); %Number = % Thresholds for rejections; Arg1 = Brain; Arg2 = Muscle; Arg3 = Eye; Arg4 = Heart; Arg5 = Line Noise; Arg6 = Channel Noise; Arg7 = Other

        EEGPreprocessTAB.num_ICA_rem = sum(EEG.reject.gcompreject);              % Overall number of ICs selected for removal
        EEGPreprocessTAB.num_musc_ICA_rem = sum(EEG.etc.ic_classification.ICLabel.classifications(:,2) > 0.9); % Number of muscle ICs selected for removal
        EEGPreprocessTAB.num_eye_ICA_rem = sum(EEG.etc.ic_classification.ICLabel.classifications(:,3) > 0.9); % Number of eye ICs selected for removal
        EEGPreprocessTAB.num_line_ICA_rem = sum(EEG.etc.ic_classification.ICLabel.classifications(:,5) > 0.9); % Number of line IC selected for removal
        EEGPreprocessTAB.num_chan_ICA_rem = sum(EEG.etc.ic_classification.ICLabel.classifications(:,6) > 0.9); % Number of channel ICs selected for removal
        
        % Step 7.3: Remove selected components
        EEG = pop_subcomp(EEG, [], 0);  % Remove components flagged by ICLabel
        EEG = eeg_checkset(EEG); % checking dataset integrity

    % Step 8: Interpolate removed bad channels
    EEG = pop_interp(EEG, EEG.chaninfo.removedchans(3:end), 'spherical');       % Interpolate all channels removed by FASTER
    EEG = eeg_checkset(EEG); % checking dataset integrity

%% Write .set  file for bandpower calculation
EEG = pop_saveset(EEG, 'filename', strcat(rawData(subjID).folder,'\',rawData(subjID).name(1:end-4),'_processed.set'), ...
      'check', 'on', 'savemode', 'onefile', 'version', '7.3'); 

FilesProcessed(subjID,:) = EEG.setname; % Log the file that was just processed
writetable(EEGPreprocessTAB, strcat(rawData(subjID).folder,'\',rawData(subjID).name(1:end-4),'preprocessing_info.csv')) % Write pre-processing info to file
end 

writematrix(FilesProcessed, strcat(cd,'\FilesProcessed_', date,'.txt'))