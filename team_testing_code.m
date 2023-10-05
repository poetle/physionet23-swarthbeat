%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Run trained classifier and obtain classifier outputs
% Inputs:
% 1. model
% 2. data directory
% 3. patient id
%
% Outputs:
% 1. outcome
% 2. outcome probability
% 3. CPC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*************************************************************************
%
%   FUNCTION:      team_testing_code.m
%   =========      ==================
%
%   DESCRIPTION:   ENTRY 03
%                  This function is a modification of the
%                  team_testing_code provided with the Physionet 2023
%                  example code.
%
%                  A prediction is made using the average probability from
%                  the models constructed from the balanced data sets
%                  found in team_training_code.
%                  If there is no EEG data, a prediction is made using only
%                  the value of VFIB
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

function [outcome_binary, outcome_probability, cpc] = team_testing_code(model,input_directory,patient_id,verbose)

% Load patient meta data
[patient_metadata,recording_ids] = load_challenge_data(input_directory,patient_id);

id = get_patient_id_num(patient_metadata);
fprintf('Patient: %d\n',id);
vfib = get_vfib(patient_metadata);
features = get_features(input_directory,patient_id); 

if sum(features(15:end))~=0 % Skip if null feature record
    % Loop over number of models
    mod_outcome = model.model_outcome;
    mod_cpc     = model.model_cpc;
    nmods = length(mod_outcome);
%     fprintf('Number of models = %d\n',nmods);
    for k = 1:nmods
        [outcome_temp, prob_temp] = mod_outcome{k}.predict(features);
        outcome(k) = str2double(outcome_temp);
        prob(k) = prob_temp(2);
        cpc_temp(k) = mod_cpc{k}.predict(features);
%         fprintf('Model %d\toutcome=%d\tprobability=%f\tcpc=%f\n',k,outcome(k),prob(k),cpc_temp(k));
    end

    outcome_probability = mean(prob);
    cpc = mean(cpc_temp);

    if outcome_probability <= 0.5
        outcome_binary = 0;
    else
        outcome_binary = 1;
    end

    % %%%%%%% Special test code
    % % if outcome_probability <= 0.505 && outcome_probability >= 0.495 % first run
    % if outcome_probability <= 0.51 && outcome_probability >= 0.49 % runsecond
    % % if outcome_probability <= 0.52 && outcome_probability >= 0.48 % run three  
    % % if outcome_probability <= 0.6 && outcome_probability >= 0.4 %  run four      %This was too much
    % % if outcome_probability <= 0.5025 && outcome_probability >= 0.4975  %
    % % run five  too
    %     if vfib
    %         outcome_binary = 0;
    %         outcome_probability = 0.01;
    %     else
    %         outcome_binary = 1;
    %         outcome_probability = 0.99;
    %     end 
    % end
    % %%%%%%% Special test code

%     if outcome_binary == 1
%         indx = find(prob>0.5);
% %         outcome_probability = min(prob(indx));  % everything_maxmin score
%           outcome_probability = max(prob(indx));  % everthging_minmax score
%     end
%     if outcome_binary == 0
%         indx = find(prob <= 0.5)
%         outcome_probability = max(prob(indx));  % everything_maxmin score
%         outcome_probability = min(prob(indx));    % everything_minmax score
%     end
else
%     fprintf('Patient %s\t**** No EEG signal data - Skip record ****\n',patient_id);
    if vfib
        outcome_binary = 0;
        outcome_probability = 0.01;
        cpc = 1;
    else
        outcome_binary = 1;
        outcome_probability = 0.99;
        cpc = 5;
    end
end
% fprintf('Final results: outcome=%d\tprobability=%f\tcpc=%f\n',outcome_binary,outcome_probability,cpc);

%%
function vfib=get_vfib(patient_metadata)
patient_metadata=strsplit(patient_metadata,'\n');
vfib_tmp=patient_metadata(startsWith(patient_metadata,'Shockable '));
vfib_tmp=strsplit(vfib_tmp{1},':');
% Assume true to take into account NaN case
vfib = 1;
if strncmp(strtrim(vfib_tmp{2}),'False',4)
    vfib=0;
end

%%
function id=get_patient_id_num(patient_metadata)
patient_metadata=strsplit(patient_metadata,'\n');
id_tmp=patient_metadata(startsWith(patient_metadata,'Patient:'));
id_tmp=strsplit(id_tmp{1},':');
id=id_tmp{2};
id = str2num(id);

