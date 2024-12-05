mkdir([path '\prepros'])
if ~iscell(filename)
    filename={filename};
end

% for each file selected, run this pipeline
for i=1:size(filename,1)

    %%% LOAD THE FILES
    EEG = loadcurry([path filename{i}], 'CurryLocations', 'True'); %load file with channel locations
    EEG = eeg_checkset(EEG); %checks the consistency of the datafields

    %%% STEP 1: REMOVE EXTRA CHANNELS AND RESAMPLE TO 1000HZ %%%
    EEG = pop_select(EEG,'nochannel',{'IO1';'IO2';'VEOG';'HEOG';'TRIGGER';'EMG';'EKG'}'); % removing unneccesary channels
    EEG.preprosTAB=table({EEG.filename},'VariableNames',{'Participant'}); % Creating report table - saving preprocessing information
    EEG.preprosTAB.Original_srate=(EEG.srate); %save the original sampling rate
    if EEG.srate>1001; EEG = pop_resample(EEG,1000); end %resample data to 1000Hz if needed

    %%% STEP 2: FILTER %%%
    EEG.data=double(EEG.data); 
    low=0.1; high=80; 
    % Setting parameters to construct the filter - 2nd order butterworth
    [x_filthigh,y_filthigh]=butter(2,low./(EEG.srate/2),'high'); %Highpass
    [x_filtlow,y_filtlow]=butter(2,high./(EEG.srate/2)); %Lowpass
    [x_filtstop,y_filtstop]=butter(4,[59./(EEG.srate/2) 61./(EEG.srate/2)],'stop'); %line noise
    for elec=1:size(EEG.data,1)
        EEG.data(elec,:)=filtfilt(x_filthigh,y_filthigh,double(EEG.data(elec,:)));
        EEG.data(elec,:)=filtfilt(x_filtlow,y_filtlow,double(EEG.data(elec,:)));
        EEG.data(elec,:)=filtfilt(x_filtstop,y_filtstop,double(EEG.data(elec,:)));
    end
    EEG.data=double(EEG.data);
    EEG.preprosTAB.filter=[low,high];
    %%% STEP 3: REMOVE BAD CHANNELS %%%
    Chan64=EEG.chanlocs; %channel location - save montage
    eyechan=find(~strcmpi('M1',{EEG.chanlocs.labels}) & ~strcmpi('M2',{EEG.chanlocs.labels}) & ~strcmpi('FP1',{EEG.chanlocs.labels})...
        & ~strcmpi('FP2',{EEG.chanlocs.labels}) & ~strcmpi('FPz',{EEG.chanlocs.labels})) ; %disregards most frontal channels - eye blinks,and M1 and M2 channels
    % Marking channels for removal with 1-6 Hz activation bigger than 2.5 zscores
    [EEG, rejelecs] = pop_rejchan(EEG,'elec',eyechan,'threshold',2.5,'measure','spec','norm','on','freqrange',[1 6]);
    EEG.rejelecs={Chan64([rejelecs]).labels}; %save rejected electrodes
    % Repeat for line noise (60Hz)
    Chan=EEG.chanlocs;
    eyechan=find(~strcmpi('M1',{EEG.chanlocs.labels}) & ~strcmpi('M2',{EEG.chanlocs.labels}) & ~strcmpi('FP1',{EEG.chanlocs.labels})...
        & ~strcmpi('FP2',{EEG.chanlocs.labels}) & ~strcmpi('FPz',{EEG.chanlocs.labels})) ; %disregards most frontal channels - eye blinks and M1 and M2 channels

    [EEG, rejelecs] = pop_rejchan(EEG,'elec',eyechan,'threshold',2.5,'measure','spec','norm','on','freqrange',[58 62]);
    EEG.rejelecs=[EEG.rejelecs {Chan([rejelecs]).labels}];

    % report rejelecs
    EEG.preprosTAB.RejElecs={string([EEG.rejelecs])}; %save the rejected electrodes as cell
    EEG.preprosTAB.RejElecsSum={size(EEG.rejelecs,2)}; %sum = how many electrodes were removed
    %%% STEP 4: ICA %%%
    EEG.data=double(EEG.data);
    EEG = pop_runica(EEG,'icatype','fastica','approach','symm','g','tanh','stabilization','on'); % running fast-ICA (dependes of FAST ICA EEGLAB addon)
    EEG = pop_iclabel(EEG, 'default'); % automatic ICA components labeling - by EEGLAB
    % Rejecting components based on ICA
    if sum(EEG.etc.ic_classification.ICLabel.classifications(:,3)>0.9,'all')>0 % IC #0-3, if 90% probability/tolerance for strength of covariance
        RejComp=find(EEG.etc.ic_classification.ICLabel.classifications(:,3)>0.9); % counting bad eye components for reporting
        TotalRejComp=sum(EEG.etc.ic_classification.ICLabel.classifications(:,3)>0.9); % counting bad eye components for reporting
        EEG = pop_subcomp(EEG,[find(any(EEG.etc.ic_classification.ICLabel.classifications(:,3)>0.9,2))] , 0); % removing eye components idendified with certinty above 90%
        EEG.preprosTAB=[EEG.preprosTAB cell2table({RejComp'},'VariableNames',{[strrep(EEG.etc.ic_classification.ICLabel.classes{3},' ','_') '_comps_removed']})]; % saving
        EEG.preprosTAB.Total_eyes_comps_removed=[sum(TotalRejComp)];
    else
        EEG.preprosTAB.Eye_comps_removed= 0;
        EEG.preprosTAB.Total_eyes_comps_removed=0;
    end
    % Reporting potentially bad ICA comps  - line noise components
    if sum(EEG.etc.ic_classification.ICLabel.classifications(:,5)>0.9,'all')>0
        EEG.preprosTAB=[EEG.preprosTAB cell2table(mat2cell(find(EEG.etc.ic_classification.ICLabel.classifications(:,5)>0.9),1),...
            'VariableNames',{[strrep(EEG.etc.ic_classification.ICLabel.classes{5},' ','_') '_suspected_comps']})];
    else
        EEG.preprosTAB.Line_Noise_suspected_comps= 0;
    end
    % Reporting suspected channel noise components
    if sum(EEG.etc.ic_classification.ICLabel.classifications(:,6)>0.9,'all')>0
        EEG.preprosTAB=[EEG.preprosTAB cell2table(mat2cell(find(EEG.etc.ic_classification.ICLabel.classifications(:,6)>0.9),1),...
            'VariableNames',{[strrep(EEG.etc.ic_classification.ICLabel.classes{6},' ','_') '_suspected_comps']})];
    else
        EEG.preprosTAB.Channel_Noise_suspected_comps= 0; 
    end
    EEG= pop_interp(EEG, Chan64, 'spherical'); % intepolating removed and missing channels
    EEG = eeg_checkset(EEG); % checking dataset integrity
    EEG.filename=[EEGA.filename(1:end-4) '.set']; EEG.setname=EEG.filename(1:end-4);
    EEG.filepath=[];
    EEG = pop_saveset(EEG, 'filename',EEG.filename,'filepath',strrep([path 'prepros\'],'\','\\'),'check', 'on','savemode','onefile','version','7.3'); % saving preprocessed subject datafile
end
