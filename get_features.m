function features=get_features(input_directory,patient_id)
%*************************************************************************
%
%   FUNCTION:      get_features.m
%   =========      ==================
%
%   DESCRIPTION:   ENTRY 03
%                  This function is a modification of the function provided
%                  by the Physionet 2023 example code.
%
%                  Starting from the latest EEG, it examines the quality of
%                  the EEG data to determine the "best" last EEG data.
%                  It only uses 5 minutes of data extracted from the middle
%                  of the hour recording.
%                  The data is resampled to 100 Hz.
%                  Bipolar signals are calculated by subtracting adjacent
%                  unipolar signals.
%                  Patient attributes for age and VFIB are used in
%                  conjunction with the EEG features returned for
%                  get_eeg_features.
%                  There was not time to incoporate additional information
%                  from the ECG data.
%                
%
%   COPYWRITE:     Allan R. Moser, Lys Kang, Jackie Le
%   ==========     Swarthmore College
%                  Engineering Department
%                  Swarthmore, PA  19081
%
%   DATE CREATED:  08-20-2023
%   =============
%
%   LAST CHANGED:  08-21-2023
%   =============
%
%**************************************************************************

% tstart = tic;   % Start timer

[patient_metadata,recording_ids] = load_challenge_data(input_directory,patient_id);
num_recordings = length(recording_ids);

% Extract patient features
id       = get_patient_id_num(patient_metadata);
age      = get_age(patient_metadata);
sex      = get_sex(patient_metadata);
rosc     = get_rosc(patient_metadata);
ohca     = get_ohca(patient_metadata);
vfib     = get_vfib(patient_metadata);
ttm      = get_ttm(patient_metadata);
hospital = get_hospital(patient_metadata);  % Add code for get_hospital

% Use one-hot encoding for sex
if strncmp(sex,'Male',4)
    male=1;
    female=0;
    other=0;
elseif strncmp(sex,'Female',4)
    male=0;
    female=1;
    other=0;
else
    male=0;
    female=0;
    other=1;
end

% patient_features=[hospital age male female other rosc ohca vfib ttm];
patient_features=[hospital age male female other rosc ohca vfib ttm];

%Extract EEG features
channels = {'Fp1','Fp2','F3','F4','F7','F8','Fz','C3','C4','Cz', ...
            'T3','T4','T5','T6','P3','P4','Pz','O1','O2'};
group='EEG';

if num_recordings>0
    % Check to see if there is any eeg data
    for j = 1:num_recordings
       recording_id = recording_ids{j};
       recording_location = fullfile(input_directory,patient_id,sprintf('%s_%s',recording_id,group));
       if exist([recording_location '.hea'],'file')>0 & exist([recording_location '.mat'],'file')>0
           break
       end
       % If it doesn't break out of this, there is no EEG data
       fprintf('No EEG recordings\n');
       features_eeg = zeros(1,168);
       hour_first = 0;
       hour_last = 0;
       hour_del = hour_last - hour_first;
       patient_hours = [hour_first hour_last hour_del];
       features_ecg = [];
       features = [age vfib features_eeg];
       return
    end

    tmp = strsplit(recording_ids{1},'_');
    hour_first = str2num(tmp{3});
    old_hour = -1;
    for j = num_recordings:-1:1
       recording_id = recording_ids{j};
       tmp = strsplit(recording_id,'_');
       hour = str2num(tmp{3});
       if (hour == old_hour)    % Skipping over duplicate hours
           old_hour = hour;
%            fprintf('Skipping:%d\t%s\t%d\n',j,recording_id,hour);
           continue;
       end
       if hour > 72
           old_hour = hour;
%            fprintf('Greater than hour 72:%d\t%s\t%d\n',j,recording_id,hour);                              % Comment out print statement
           continue;
       end
       old_hour = hour;
%        fprintf('j=%d\trecording_id=%s\t,hour=%d\n',j,recording_id,hour);                                  % Comment out print statement
       recording_location = fullfile(input_directory,patient_id,sprintf('%s_%s',recording_id,group));
       if exist([recording_location '.hea'],'file')>0 & exist([recording_location '.mat'],'file')>0
%            fprintf('recording = %d\trecording_location: %s\n',j,recording_location);                      % Comment out print statement
           [signal_data,sampling_frequency,signal_channels]=load_recording(recording_location);
           
           utility_frequency=get_utility_frequency(recording_location);
           [signal_data, ~] = reduce_channels(signal_data, channels, signal_channels);
           [signal_data, sampling_frequency] = preprocess_data(signal_data, sampling_frequency, utility_frequency);

           % Check for no data after trying to extract 5 minute segment
           no_data = isempty(signal_data);
           if no_data
%                fprintf('Not enough data: j=%d\trecording_id=%s\t,hour=%d\n',j,recording_id,hour);            % Comment out print statement
               sig_qual = 0;
               continue;
           end

           [sig_qual,nc,num_zeros,sdev] = check_quality(signal_data); 
           if (~sig_qual)
                 % Comment out print statement
%                fprintf('Failed quality test:j=%d\trecording_id=%s\t,hour=%d\tchannel=%d\tnum_zeros=%d\tsdev=%10.2e\n',...
%                         j,recording_id,hour,nc,num_zeros,sdev);
               continue;
           end

           % Check to see if bad data
           if ~isempty(find(isnan(signal_data)))
%                fprintf('Values in signal_data that are NaN: %d\t%s\t%d\n',j,recording_id,hour);   % Comment out print statement
               sig_qual = 0;
               continue
           end
       else
%            fprintf('Header or mat file do not exist for this hour: %d\t%s\t%d\n',j,recording_id,hour);
           continue;
       end
       break;
   end
   if (~sig_qual || no_data)
        if ~sig_qual
            fprintf('All Failed quality test: j=%d\trecording_id=%s\t,hour=%d\n',j,recording_id,hour);      % Comment out print statement
        end
        if no_data
%             fprintf('Not enough data: j=%d\trecording_id=%s\t,hour=%d\n',j,recording_id,hour);              % Comment out print statement
        end
        hour_last = hour;
        features_eeg =  zeros(1,168);
   else
       hour_last = hour;
%        fprintf('Last record read is %d, sampling frequency = %d, number of samples = %d\n',j,sampling_frequency,length(signal_data));   % Comment out print statement
       % Convert to bipolar montage: F3-P3 and F4-P4
       data( 1,:)=signal_data(05,:)-signal_data(01,:);   %  1: F7 - Fp1
       data( 2,:)=signal_data(11,:)-signal_data(05,:);   %  2: T3 - F7
       data( 3,:)=signal_data(13,:)-signal_data(11,:);   %  3: T5 - T3
       data( 4,:)=signal_data(18,:)-signal_data(13,:);   %  4: O1 - T5
       data( 5,:)=signal_data(06,:)-signal_data(02,:);   %  5: F8 - Fp2
       data( 6,:)=signal_data(12,:)-signal_data(06,:);   %  6: T4 - F8
       data( 7,:)=signal_data(14,:)-signal_data(12,:);   %  7: T6 - T4
       data( 8,:)=signal_data(19,:)-signal_data(14,:);   %  8: O2 - T6
       data( 9,:)=signal_data(03,:)-signal_data(01,:);   %  9: F3 - Fp1
       data(10,:)=signal_data(08,:)-signal_data(03,:);   % 10: C3 - F3
       data(11,:)=signal_data(15,:)-signal_data(08,:);   % 11: P3 - C3
       data(12,:)=signal_data(18,:)-signal_data(15,:);   % 12: O1 - P3
       data(13,:)=signal_data(04,:)-signal_data(02,:);   % 13: F4 - Fp2
       data(14,:)=signal_data(09,:)-signal_data(04,:);   % 14: C4 - F4
       data(15,:)=signal_data(16,:)-signal_data(09,:);   % 15: P4 - C4
       data(16,:)=signal_data(19,:)-signal_data(16,:);   % 16: O2 - P4
       data(17,:)=signal_data(10,:)-signal_data(07,:);   % 17: Cz - Fz
       data(18,:)=signal_data(17,:)-signal_data(10,:);   % 18: Pz - Cz
       dat = data';
       % Finally, get eeg features!
       features_eeg = get_eeg_features(dat,sampling_frequency);
   end
else
    fprintf('No recordings\n');
    features_eeg = zeros(1,168);
    hour_first = 0;
    hour_last = 0;
end
   

% Combine the features from the patient metadata and the recording data and metadata.

hour_del = hour_last - hour_first;
patient_hours = [hour_first hour_last hour_del];
features_ecg = [];
features = [age vfib features_eeg];
% toctime = toc(tstart);     %stop timer
% fprintf('Elaspsed time = %f\n',toctime)

%%
function id=get_patient_id_num(patient_metadata)
patient_metadata=strsplit(patient_metadata,'\n');
id_tmp=patient_metadata(startsWith(patient_metadata,'Patient:'));
id_tmp=strsplit(id_tmp{1},':');
id=id_tmp{2};
id = str2num(id);


function hospital = get_hospital(patient_metadata)
patient_metadata=strsplit(patient_metadata,'\n');
hospital_tmp=patient_metadata(startsWith(patient_metadata,'Hospital:','IgnoreCase',true));
hospital_tmp=strsplit(hospital_tmp{1},':');
hospital=strtrim(hospital_tmp{2});
hospital = double(upper(hospital)) - double('A') +1; % Convert letters to numbers

function age=get_age(patient_metadata)
patient_metadata=strsplit(patient_metadata,'\n');
age_tmp=patient_metadata(startsWith(patient_metadata,'Age:'));
age_tmp=strsplit(age_tmp{1},':');
age=str2double(age_tmp{2});

function sex=get_sex(patient_metadata)
patient_metadata=strsplit(patient_metadata,'\n');
sex_tmp=patient_metadata(startsWith(patient_metadata,'Sex:'));
sex_tmp=strsplit(sex_tmp{1},':');
sex=strtrim(sex_tmp{2});

function rosc=get_rosc(patient_metadata)
patient_metadata=strsplit(patient_metadata,'\n');
rosc_tmp=patient_metadata(startsWith(patient_metadata,'ROSC:'));
rosc_tmp=strsplit(rosc_tmp{1},':');
rosc=str2double(rosc_tmp{2});

function ohca=get_ohca(patient_metadata)
patient_metadata=strsplit(patient_metadata,'\n');
ohca_tmp=patient_metadata(startsWith(patient_metadata,'OHCA:'));
ohca_tmp=strsplit(ohca_tmp{1},':');
if strncmp(strtrim(ohca_tmp{2}),'True',4)
    ohca=1;
elseif strncmp(strtrim(ohca_tmp{2}),'False',4)
    ohca=0;
else
    ohca=nan;
end

function vfib=get_vfib(patient_metadata)
patient_metadata=strsplit(patient_metadata,'\n');
vfib_tmp=patient_metadata(startsWith(patient_metadata,'Shockable '));
vfib_tmp=strsplit(vfib_tmp{1},':');
if strncmp(strtrim(vfib_tmp{2}),'True',4)
    vfib=1;
elseif strncmp(strtrim(vfib_tmp{2}),'False',4)
    vfib=0;
else
    vfib=nan;
end

function ttm=get_ttm(patient_metadata)
patient_metadata=strsplit(patient_metadata,'\n');
ttm_tmp=patient_metadata(startsWith(patient_metadata,'TTM:'));
ttm_tmp=strsplit(ttm_tmp{1},':');
ttm=str2double(ttm_tmp{2});

function [TF,nc,num_zeros,sdev] = check_quality(data)   % This is ARM quality test
[nrows,ncols] = size(data);
nchannels = min(nrows,ncols);
nsamps = max(nrows,ncols);
max_zeros = nsamps/2;
min_sd = 0.001;
num_zeros = NaN;
sdev = NaN;
TF = 1;  % Good quality
for nc = 1:nchannels
    num_zeros = length(find(data(nc,:) == 0));
    if num_zeros > max_zeros
        TF = 0;
        break
    end
    sdev = std(data(nc,:));
    if sdev < min_sd
        TF = 0;
        break
    end
end

function utility_frequency=get_utility_frequency(recording_info)
header_file=strsplit([recording_info '.hea'],'/');
header_file=header_file{end};
header=strsplit(fileread([recording_info '.hea']),'\n');
utility_tmp=header(startsWith(header,'#Utility'));
utility_tmp=strsplit(utility_tmp{1},':');
utility_frequency=str2double(utility_tmp{2});

function [data, channels] = reduce_channels(data, channels, signal_channels)
channel_order=signal_channels(ismember(signal_channels, channels));
data=data(ismember(signal_channels, channels),:);
data=reorder_recording_channels(data, channel_order, channels);

function reordered_signal_data=reorder_recording_channels(signal_data, current_channels, reordered_channels)
if length(current_channels)<length(reordered_channels)
    for i=length(current_channels)+1:length(reordered_channels)
        current_channels{i}='';
    end
end
if sum(cellfun(@strcmp, reordered_channels, current_channels))~=length(current_channels)
    indices=[];
    for j=1:length(reordered_channels)
        if sum(strcmp(reordered_channels{j},current_channels))>0
            indices=[indices find(strcmp(reordered_channels{j},current_channels))];
        else
            indices=[indices nan];
        end
    end
    num_channels=length(reordered_channels);
    num_samples=size(signal_data,2);
    reordered_signal_data=zeros(num_channels,num_samples);
    for j=1:num_channels
        if ~isnan(indices(j))
            reordered_signal_data(j,:)=signal_data(indices(j),:);
        end
    end
else
    reordered_signal_data=signal_data;
end

function [rescaled_data,sampling_frequency,channels]=load_recording(recording_location)
header_file=strsplit([recording_location '.hea'],'/');
header_file=header_file{end};
header=strsplit(fileread([recording_location '.hea']),'\n');
header(cellfun(@(x) isempty(x),header))=[];
header(startsWith(header,'#'))=[];
recordings_info=strsplit(header{1},' ');
record_name=recordings_info{1};
num_signals=str2double(recordings_info{2});
sampling_frequency=str2double(recordings_info{3});
num_samples=str2double(recordings_info{4});
signal_file=cell(1,length(header)-1);
gain=zeros(1,length(header)-1);
offset=zeros(1,length(header)-1);
initial_value=zeros(1,length(header)-1);
checksum=zeros(1,length(header)-1);
channels=cell(1,length(header)-1);
for j=2:length(header)
    header_tmp=strsplit(header{j},' ');
    signal_file{j-1}=header_tmp{1};
    gain(j-1)=str2double(header_tmp{3});
    offset(j-1)=str2double(header_tmp{5});
    initial_value(j-1)=str2double(header_tmp{6});
    checksum(j-1)=str2double(header_tmp{7});
    channels{j-1}=header_tmp{9};
end
if ~length(unique(signal_file))==1
    error('A single signal file was expected for %s',header_file)
end

% Load the signal file
[val, name, T] = load_mat_v4([recording_location '.mat']);

num_channels=length(channels);
if num_channels~=size(val,1) || num_samples~=size(val,2)
    error('The header file %s is inconsistent with the dimensions of the signal file',header_file)
end

for j=1:num_channels
    if val(j,1)~=initial_value(j)
        error('The initial value in header file %s is inconsistent with the initial value for the channel',header_file)
    end
    
    if sum(val(j,:))~=checksum(j)
        error('The checksum in header file %s is inconsistent with the initial value for the channel',header_file)
    end
end

rescaled_data=zeros(num_channels,num_samples);
for j=1:num_channels
    rescaled_data(j,:)=(val(j,:)-offset(j))/gain(j);
end

%%
function [data, resampling_frequency]=preprocess_data(data, sampling_frequency, utility_frequency)
% Re-wrote most of this to resample correctly

% Define the bandpass frequencies.
passband = [0.1, 30.0];
resampling_frequency = 100;

% Resample to resamp_frequency rate
% This applies an FIR anti-aliasing filter (by default of order 10)
% Note that array has sample number first, channel number second
resamp_data = resample(data',resampling_frequency,sampling_frequency);
% Cut off the first and last bits data to avoid end effects of FIR filter
if length(resamp_data) > 50
    resamp_data = resamp_data(21:end-20,:);
end

% Demean the data
resamp_data = resamp_data - mean(resamp_data);

% Clip out middle 5 minutes if enough data
newsamp = 30000;
nsamp = length(resamp_data);
if nsamp > newsamp
    start = floor((nsamp - newsamp)/2) + 1;
    resamp_data = resamp_data(start:start+newsamp-1,:);
else nsamp < newsamp;
    resamp_data = [];
end 

data = resamp_data';

%%
function features_ecg=get_ecg_features(data, sampling_frequency)

features_ecg=[mean(data') std(data')];
