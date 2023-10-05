function model = team_training_code(input_directory,output_directory, verbose) % train_EEG_classifier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Train EEG classifiers and obtain the models
% Inputs:
% 1. input_directory
% 2. output_directory
%
% Outputs:
% 1. model: trained model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*************************************************************************
%
%   FUNCTION:      team_training_code.m
%   =========      ==================
%
%   DESCRIPTION:   ENTRY 03
%                  This function is a modification of the
%                  team_training_code provided with the Physionet 2023
%                  example code.
%
%                  An important difference is that it randomly subsamples the
%                  majority class (poor outcome) instances to obtain
%                  data sets that have equal numbers of good and poor
%                  outcome instances.  This is done multiple times. A new
%                  model is constructed for each balanced data set.
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

if verbose>=1
    disp('Finding challenge data...')
end

% Find the folders
patient_ids=dir(input_directory);
patient_ids=patient_ids([patient_ids.isdir]==1);
patient_ids(1:2)=[]; % Remove "./" and "../" paths 
patient_ids={patient_ids.name};
num_patients = length(patient_ids);

% Create a folder for the model if it doesn't exist
if ~isdir(output_directory)
    mkdir(output_directory)
end

fprintf('Loading data for %d patients...\n', num_patients)

% Extract fatures and labels

features=[];
outcomes=[];                                                               % arm change
cpcs=[];

k = 0;                                                                     % arm change
for j=1:num_patients
    if verbose>1
        fprintf('%d/%d \n',j,num_patients)
    end

    % Extract features
    patient_id=patient_ids{j};
    current_features=get_features(input_directory,patient_id);              

    if sum(current_features(13:end))~=0 % Skip if null feature record      
        k = k+1;
        features(k,:)=current_features;

        % Load data
        patient_metadata=load_challenge_data(input_directory,patient_id);
    
        % Extract labels
        current_outcome=get_outcome(patient_metadata);
        outcomes(k)=current_outcome;
        current_cpc=get_cpc(patient_metadata);
        cpcs(k)=current_cpc;
        patient_ids_keep{k} = patient_id;
        fprintf('Patient: %s\tOutcome: %d\tCPC: %d\n',patient_ids_keep{k},outcomes(k),cpcs(k))
    else 
        fprintf('**** No EEG signal data - Skip record %d ****\n',j)
    end   
end

%% Create several balanced data sets
good_indx = find(outcomes==0);
poor_indx = find(outcomes==1);
num_good = length(good_indx);
num_poor = length(poor_indx);

% Determine which class has the fewest instances
ntot_per_class = min(num_good,num_poor);
fprintf('Number good = %d\tNumber poor = %d\tNumber to use = %d\n',num_good,num_poor,ntot_per_class);
% Use same number of instances per class

% nsets = 5;  % Number of balanced datasets
nsets = 10;
for k = 1:nsets
    rand_poor_indx = randperm(num_poor);
    rand_good_indx = randperm(num_good);
    sel_poor_indx = poor_indx(rand_poor_indx(1:ntot_per_class));
    sel_good_indx = good_indx(rand_good_indx(1:ntot_per_class));
    sel_indx = [sel_good_indx,sel_poor_indx];
    sel_indx = sort(sel_indx);
    sel_features = features(sel_indx,:);
    %% 
    sel_outcomes = outcomes(sel_indx);
    sel_cpcs     = cpcs(sel_indx); 
    fprintf('Training model %d\n',k);
    % model_outcome{k} = TreeBagger(5000,sel_features,sel_outcomes);
    template = templateTree('MaxNumSplits', 1); 
    costMatrix = [0 3; 3 0]; 
    model_outcome{k} = fitcensemble(sel_features,sel_outcomes, 'Method', 'AdaBoostM1', 'Learners', template, 'Cost', costMatrix);
    model_cpc{k}     = TreeBagger(5000,sel_features,sel_cpcs,'method','regression');
end
save_model(model_outcome,model_cpc,output_directory);
disp('Done.');
end
%%

function save_model(model_outcome,model_cpc,output_directory) %save_model
% Save results.
filename = fullfile(output_directory,'model.mat');
save(filename,'model_outcome','model_cpc','-v7.3');

disp('Done.')
end

function outcome=get_outcome(patient_metadata)

patient_metadata=strsplit(patient_metadata,'\n');
outcome_tmp=patient_metadata(startsWith(patient_metadata,'Outcome:'));
outcome_tmp=strsplit(outcome_tmp{1},':');

if strncmp(strtrim(outcome_tmp{2}),'Good',4)
    outcome=0;
elseif strncmp(strtrim(outcome_tmp{2}),'Poor',4)
    outcome=1;
else
    keyboard
end

end

function cpc=get_cpc(patient_metadata)

patient_metadata=strsplit(patient_metadata,'\n');
cpc_tmp=patient_metadata(startsWith(patient_metadata,'CPC:'));
cpc_tmp=strsplit(cpc_tmp{1},':');
cpc=str2double(cpc_tmp{2});

end