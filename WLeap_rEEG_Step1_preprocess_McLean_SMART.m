 %%% WellcomeLEAP - Resting State EEG data pre-processing %%%

% Last updated: 2024/01/18
% Author: Ty Lees

% This script requires:
% EEGLAB 2023.1 or later (for the re-reference step)
% The following extensions/plugins (beyond the default)
% ERPLAB, FASTER, ICLabel, clean_rawdata, bva-io

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

cd('C:/Users/tlees/Partners HealthCare Dropbox/Ty Lees/WL_SMART_Trial/EEG Data Analyses/Baseline/eeg');     %Set your CD to where you keep the data; this will work with cloud directories if you have the desktop apps

%% Load raw files and concatenate together
rawData = dir('*_EO*_S1.vhdr');    % Search for all files for one Px ID
[~,index] = sortrows({rawData.name}.'); rawData = rawData(index); clear index

% Load all eye's open files for Px
for i = 1:4:length(rawData)   
    for j = 0:3
        % Load files in sets of 4; checking their names align
        % This IF is super clunky, but it works to check the IDs of each file
        if (rawData(i).name(6:8) == rawData(i+1).name(6:8)) & (rawData(i).name(6:8) == rawData(i+2).name(6:8)) & (rawData(i).name(6:8) == rawData(i+3).name(6:8)) & ...
                (rawData(i+1).name(6:8) == rawData(i+2).name(6:8)) & (rawData(i+1).name(6:8) == rawData(i+3).name(6:8)) & ...
                (rawData(i+2).name(6:8) == rawData(i+3).name(6:8))
            
            EEG = pop_loadbv(rawData(i+j).folder, rawData(i+j).name);
            EEG.setname = rawData(i+j).name;
            [ALLEEG, EEG] = eeg_store(ALLEEG, EEG);
        else
            error('One or more of the files that the script attempted to loaded did not have the correct S_ID')
        end
    end 

    % Append the 4 loaded files together
    ALLEEG(end+1) = pop_mergeset(ALLEEG, [1 2 3 4], 0);
    ALLEEG(end).setname = strcat(rawData(i).name(1:end-5), '_merged');
        
    % Save merged dataset as .set file
       %if not(isfolder('Merged Data'))
       %     mkdir 'Merged Data'
       %end 
    
    pop_saveset(ALLEEG(end), 'filename', strcat(rawData(i).folder,'\EO Merged Data\',rawData(i).name(1:end-9),'_S1_merged.set'), ...
                'check', 'on', 'savemode', 'onefile', 'version', '7.3');
    
    ALLEEG = []; EEG=[]; CURRENTSET=[]; % Clear all loaded files before next loop
end     

%% Load file and pre-process data
eeglab

mergedData = dir('*\*_merged.set');
FilesProcessed = strings(length(mergedData),1);

% Load Stanford file for virtual channel locations
Stanford = pop_loadset('Stanford_locs.set',...
                   'C:\Users\tlees\Partners HealthCare DropBox\Ty Lees\Postdoc - McLean\Research Projects\WellcomeLEAP_personal\Test Data for Pipelines\McLean Test Data\');

for subjID = 1:length(mergedData) 

    % Step 1: Import data; function depends on data format
    loadName = strcat(mergedData(subjID).folder,'\',mergedData(subjID).name);
    EEG = pop_loadset(loadName);  EEG.setname = mergedData(subjID).name(1:end-4);    % Load and name merged file
    EEGPreprocessTAB = table({EEG.setname},'VariableNames',{'Participant'});      % Record participant ID and create table for processing data
    EEG = eeg_checkset(EEG); % Check dataset integrity

    % Step 2: Add Channel information and remove unwanted channels    
    EEG = pop_chanedit(EEG, 'append', 95, 'changefield', {96,'labels','1'}, ...         % Add ref electrode ('1') location information and set as a non-data channel
                      'changefield',{96,'sph_radius','1'}, 'changefield',{96,'sph_theta','30'},'changefield',{96,'sph_phi','90'},...
                      'convert',{'sph2all'}, 'changefield',{96,'datachan',0}); 
    EEG = pop_chanedit(EEG, 'changefield', {96, 'type', 'REF'});    % Set electrode 1 type to 'reference'
    % No additional channels are recorded and need to be removed

    % Step 3: Resample to 1000 Hz if needed
    EEGPreprocessTAB.Orig_srate = EEG.srate;    % Record original sample rate

    if EEG.srate > 1000          
        EEG = pop_resample(EEG, 1000);
    elseif EEG.srate < 1000
        warning('Sample rate below target; function is trying to up-sample')
        fid = fopen(strcat(loadName(1:end-11),'_Resample_warning.txt'), 'wt');
            fprintf(fid, 'Recording sample rate is below target, and function is trying to up-sample');
        fclose(fid);
    end

    % Step 4: Remove bad channels - FASTER method - Requires all channels to have location info
    EEG = pop_reref(EEG, '40', 'keepref' ,'on', 'refloc', struct('labels',{'1'}, 'sph_radius', {1}, 'sph_theta', {30}, 'sph_phi', {90}, ...     % Temporarily re-ref to Pz, add Cz back to data
                  'theta', {-30}, 'radius', {0}, 'X', {5.3029e-17}, 'Y', {3.0616e-17}, 'Z', {1}, ...
                  'type', {'REF'}, 'ref', {''}, 'urchan', {96}, 'datachan', {0}));

    list_properties = channel_properties(EEG, 1:EEG.nbchan, 39);     % Computes bad channel coefficents; Arguments are: Data input, Which Chans, Ref Chan ID   
    badChanindices = min_z(list_properties);
    badChans = find(badChanindices == 1);   % Identifies bad channels
    badChans = badChans(badChans ~= 39);     % Removes the reference from list so it's not excluded   					
    
    EEGPreprocessTAB.num_bad_chans = numel(badChans);       % Record number of bad channels identified
    EEGPreprocessTAB.badchan_indices = num2str(reshape(badChans',1,[]));    % Record index of each bad channel
    
    EEG = pop_select(EEG,'nochannel', badChans);    % Remove bad channels from recording   

    % Step 5: Filter data
        % Step 5.1: Bandpass filter at 0.1 - 80 Hz
        EEG  = pop_basicfilter(EEG, 1:EEG.nbchan, 'Boundary', 'boundary', 'Cutoff', [0.1 80], 'Design', 'butter', 'Filter', 'bandpass', 'Order', 2, 'RemoveDC', 'on');
        % Step 5.2: Notch filter @ 60 Hz; needs to be changed to 50 Hz for EU sites
        EEG  = pop_basicfilter(EEG, 1:EEG.nbchan, 'Boundary', 'boundary', 'Cutoff',  60, 'Design', 'notch', 'Filter', 'PMnotch', 'Order', 180); 

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

    % Step 7: Interpolate removed bad channels
    EEG = pop_interp(EEG, EEG.chaninfo.removedchans(1:end), 'spherical');
    EEG = eeg_checkset(EEG); % checking dataset integrity
    
    % Step 8: Re-reference to common average
    EEG = pop_reref(EEG, []);
    EEG = eeg_checkset(EEG); % checking dataset integrity

    % Step 9: Interpolate to UCSD electrode positions
    EEG = pop_interp(EEG, Stanford.chanlocs, 'spherical');
    EEG = pop_select(EEG, 'channel',{'Fz','F3','F7','FT9','FC5','FC1','C3','TP9','CP5','CP1','Pz','P3','P7','O1','Oz','O2', ...
                                     'P4','P8','TP10','CP6','CP2','C4','T8','FT10','FC6','FC2','F4','F8','AF7','AF3','AFz', 'F1', ...
                                     'F5','FT7','FC3','C1','C5','TP7','CP3','P1','P5','PO7','PO3','POz','PO4','PO8','P6','P2', ...
                                     'CPz','CP4','TP8','C6','C2','FC4','FT8','F6','AF8','AF4','F2','FCz','Cz','Fp1','T7','Fp2'});  %Select the virtual 10-10 channels; removes physical channels
    EEG = eeg_checkset(EEG); % checking dataset integrity

%% Write .set  file for bandpower calculation
EEG = pop_saveset(EEG, 'filename', strcat(mergedData(subjID).folder,'\',mergedData(subjID).name(1:end-11),'_processed.set'), ...
      'check', 'on', 'savemode', 'onefile', 'version', '7.3'); 

FilesProcessed(subjID,:) = EEG.setname; % Log the file that was just processed
writetable(EEGPreprocessTAB, strcat(mergedData(subjID).folder,'\',mergedData(subjID).name(1:end-11),'_preprocessing_info.csv')) % Write pre-processing info to file
end 

writematrix(FilesProcessed, strcat(cd,'\WL_SMART_FilesProcessed_', string(date),'.txt'))