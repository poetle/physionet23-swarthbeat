%*************************************************************************
%
%   FUNCTION:      get_eeg_features.m
%   =========      ==================
%
%   DESCRIPTION:   ENTRY 05
%                  This entry uses features based on signifcance determined
%                  using the estimated feature significance from the
%                  Treebagger, similar to what was done for Entry 02.
%                  The difference is that classes of features rather than
%                  individual features are used. We believe this will
%                  reduce the overtraining apparent in Entry 02.
%                  The classes of features used are:
%                  (1) Bandpower in delta, theta, alpha, and beta bands
%                  (2) The slope and r^2 found by fitting a linear function
%                      to the log of the PSD for the delta band
%                  (3) The ratio of the bandpower for the delta-to-theta
%                      and delta-to-alpha bands
%                  (4) The mean of the magnituded squared coherence
%                  estimate for the delta, theta, and alpha bands.
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
function eeg_features=get_eeg_features(dat,sampling_frequency)

    fs = sampling_frequency; % Shorter variable name for sampling frequency

    % Label individual signals to help interpretation 
    % Fp1F7 = dat(:,1);  F7T3 = dat(:,2);  T3T5 = dat(:,3);  T5O1 = dat(:,4);
    % F8T4  = dat(:,5);  F8T4 = dat(:,6);  T4T6 = dat(:,7);  T6O2 = dat(:,8);
    % Fp1F3 = dat(:,9);  F3C3 = dat(:,10); C3P3 = dat(:,11); P3O1 = dat(:,12);
    % Fp2F4 = dat(:,13); F4C4 = dat(:,14); C4P4 = dat(:,15); P4O2 = dat(:,16);
    % FzCz  = dat(:,17); CzPz = dat(:,18);

    % Get power for signal, delta, theta, alpha, beta, theta-alpha-beta for
    % all EEG lead
    bp_total = bandpower(dat,fs,[0,26]);
    bp_delta = bandpower(dat,fs,[1,3]);
    bp_theta = bandpower(dat,fs,[3,6]);
    bp_alpha = bandpower(dat,fs,[6,10]);
    bp_beta  = bandpower(dat,fs,[10,26]);

    % For opposite side of brain calculations the order is:
    % 1: Fp1 - F7 / Fp2 - F8
    % 2: F7  - T3 / F8  - T4
    % 3: T3  - T5 / T4  - T6
    % 4: T5  - O1 / T6  - O2
    % 5: Fp1 - F3 / Fp2 - F4
    % 6: F3  - C3 / F4  - C4
    % 7: C3  - P3 / C4  - P4
    % 8: P3  - O1 / P4  - O2

       
    % Find the power spectrum for all EEG leads
    [psd,pf] = pspectrum(dat,fs);
    % Find index into psd for individual bands
    pindx_delta = find(pf >  1 & pf <  3); % delta band
    pindx_theta = find(pf >  3 & pf <  6); % theta band
    pindx_alpha = find(pf >  6 & pf < 10); % alpha band
    pindx_beta  = find(pf > 10 & pf < 26); % beta band
 
    % Cross-spectral coherence for opposite sides of brain
    if length(dat) > 1024
        for i = 1:4 
            [cpsd(:,i),cf] = mscohere(dat(:,i),dat(:,i+4),1024,[],[],fs);
            k=i+4; m = i+8;
            [cpsd(:,k),cf] = mscohere(dat(:,m),dat(:,m+4),1024,[],[],fs);
        end
    else
        for i = 1:8
            cpsd(:,i) = 0;
        end
    end
    cindx_delta = find(cf >  1 & cf <  3); % delta band
    cindx_theta = find(cf >  3 & cf <  6); % theta band
    cindx_alpha = find(cf >  6 & cf < 10); % alpha band
    cindx_beta  = find(cf > 10 & cf < 26); % beta band

    % Process delta band
    freq = pf(pindx_delta);
    dbps  = 10*log10(psd(pindx_delta,:));
    for i = 1:18
        [coefs,parm] = polyfit(freq,dbps(:,i),1);
        rsq_delta(i) = 1 - parm.normr^2 / norm(dbps(:,i) - mean(dbps(:,i)))^2;
        slp_delta(i) = coefs(1);
        rat_delta_total(i)  = bp_delta(i)/bp_total(i);
        rat_delta_theta(i)  = bp_delta(i)/bp_theta(i);
        rat_delta_alpha(i)  = bp_delta(i)/bp_alpha(i);
        rat_delta_beta(i)   = bp_delta(i)/bp_beta(i);
    end
    cohere_delta = mean(cpsd(cindx_delta,:));
    
    % Process theta band
    freq = pf(pindx_theta);
    dbps  = 10*log10(psd(pindx_theta,:));
    for i = 1:18
        [coefs,parm] = polyfit(freq,dbps(:,i),1);
        rsq_theta(i) = 1 - parm.normr^2 / norm(dbps(:,i) - mean(dbps(:,i)))^2;
        slp_theta(i) = coefs(1);
        rat_theta_total(i)  = bp_theta(i)/bp_total(i);
        rat_theta_alpha(i)  = bp_theta(i)/bp_alpha(i);
        rat_theta_beta(i)   = bp_theta(i)/bp_beta(i);
    end
    cohere_theta = mean(cpsd(cindx_theta,:));
    
    % Process alpha band
    freq = pf(pindx_alpha);
    dbps  = 10*log10(psd(pindx_alpha,:));
    for i = 1:18
        [coefs,parm] = polyfit(freq,dbps(:,i),1);
        rsq_alpha(i) = 1 - parm.normr^2 / norm(dbps(:,i) - mean(dbps(:,i)))^2;
        slp_alpha(i) = coefs(1);
        rat_alpha_total(i)  = bp_alpha(i)/bp_total(i);
        rat_alpha_beta(i)   = bp_alpha(i)/bp_beta(i);
    end
    cohere_alpha = mean(cpsd(cindx_alpha,:));

    % vchan = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18];
    % xchan = [1 2 3 4 5 6 7 8];

    %t1:  F7-T3 / F8-T4	T5-O1 / T6-O2	F3-C3 / F4-C4	P3-O1 / P4-O2
    % vchan = [2, 4, 6, 8, 10, 12, 14, 16];
    % xchan = [2, 4, 6, 8];

    %t2:  F7-T3 / F8-T4	T3-T5 / T4-T6	T5-O1 / T6-O2	Fp1-F3 / FP2-F4    F3-C3 / F4-C4   C3-P3 / C4-P4
    % vchan = [2, 3, 4, 6, 7, 8, 9, 10, 11, 13, 14, 15];
    % xchan = [2, 3, 4, 5, 6, 7];

    %t3:  F7-T3 / F8-T4	T3-T5 / T4-T6	Fp1-F3 / FP2-F4	  F3-C3 / F4-C4	   C3-P3 / C4-P4
    % vchan = [2, 3, 6, 7, 9, 10, 11, 13, 14, 15];
    % xchan = [2, 3, 5, 6, 7];

    %t4:  F7-T3 / F8-T4	T3-T5 / T4-T6	Fp1-F3 / FP2-F4	   C3-P3 / C4-P4    P3-O1 / P4-O2
    % vchan = [2, 3, 6, 7, 9, 11, 12, 13, 15, 16];
    % xchan = [2, 3, 5, 7, 8];

    %t5:  F7-T3 / F8-T4	T5-O1 / T6-O2	F3-C3 / F4-C4	C3-P3 / C4-P4	P3-O1 / P4-O2
    % vchan = [2, 4, 6, 8, 10, 11, 12, 14, 15, 16];
    % xchan = [2, 4, 6, 7, 8];

    %t6:  F7-T3 / F8-T4	T3-T5 / T4-T6	C3-P3 / C4-P4
    vchan  = [2, 3, 6, 7, 11, 15];
    xchan = [2, 3, 7];

  
    eeg_features =  [bp_delta(vchan),bp_theta(vchan),bp_alpha(vchan),bp_beta(vchan), ... 
                     slp_delta(vchan),rsq_delta(vchan),rat_delta_theta(vchan),rat_delta_alpha(vchan), ...
                     cohere_delta(xchan),cohere_theta(xchan),cohere_alpha(xchan)];