
function [meanr_PR,meanr_QRS,meanr_RR,msg] = time_measures (ecg,fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Created:           Claudia Nagel (30.06.2020)
%    Last modified:     Claudia Nagel (30.06.2020)
%    Ver. 1.0.0
%    Institute of Biomedical Engineering
%    Karlsruhe Institute of Technology
%
%    http://www.ibt.kit.edu
%
%    Copyright 2000-2020 - All rights reserved.
%
% ------------------------------------------------------
% If filters an examplary 12 lead ECG (basline removal, frequency filtering
% and isoline correction). Subsequently, the annotation process is
% performed and the resulting wave boundaries for the P and T wave and the
% QRS complex are visulized with markers on the signal. 

n=10000; %number of samples i want to visualize
t = (0:length(ecg.')-1) / fs; 
figure; 
plot(t(1:n),ecg(1:n)); 
title('Unfiltered ECG Signal Lead I');
xlabel('Time (s)'); ylabel('Amplitude (mV)');
axis tight
hold all; 

%%% filtering of the 12-lead ECG
% % Remove baseline wander
% % usage: [filtered_signal,baseline]=ECG_Baseline_Removal(signal,samplerate,window_length,overlap)
[ecg_filtered_baseline,~] = ECG_Baseline_Removal(ecg,fs,1,0.5);

% visualizing waveforms
figure
subplot(3,1,1);
plot(t(1:n),ecg_filtered_baseline(1:n,1)); 
title('Baseline Removal');
xlabel('Time (s)'); ylabel('Amplitude (mV)');
axis tight

% filter noise frequencies
% frequencies are already optimized for ECG signals (literature values):
% Lowpass: 120 Hz, Highpass: 0.3 Hz, Bandstop (49-51 Hz)
[ecg_filtered_frq] = ECG_High_Low_Filter(ecg,fs,1,150);
ecg_filtered_frq=Notch_Filter(ecg_filtered_frq,fs,50,1);

% visualizing waveforms
subplot(3,1,2);
plot(t(1:n),ecg_filtered_frq(1:n,1)); 
title('Highpass, Lowpass, Bandstop');
xlabel('Time (s)'); ylabel('Amplitude (mV)');
axis tight

%%% isoline correction
% usage: [filteredsignal,offset,frequency_matrix,bins_matrix]=Isoline_Correction(signal,varargin)
[ecg_filtered_isoline,offset,~,~]=Isoline_Correction(ecg_filtered_frq);

% visualizing waveforms
subplot(3,1,3);
plot(t(1:n),ecg_filtered_isoline(1:n,1)); 
title('Isoline Correction');
xlabel('Time (s)'); ylabel('Amplitude (mV)');
axis tight

%%% Feature calculation
% produce FPT Table
% usage: [FPT_MultiChannel,FPT_Cell]=Process_ECG_Multi(signal,samplerate,varargin)
[FPT_MultiChannel,FPT_Cell]=Annotate_ECG_Multi(ecg_filtered_isoline,fs);

% extract FPTs for Channel 1 (Lead I):
FPT_LeadI = FPT_Cell{1,1};
fig=ecg_filtered_isoline(1:n,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Created:           Francesca Cavallo (25.11.2024)
%
% 
% The following part of the script returns the fiducial points of each wave that are useful to identify the time ranges of the signal
%

% R peaks index
R_peaks_i = reshape(FPT_LeadI(:, 6)', 1, []);

% visualize fiducial points
figure; 
title('Filtered ECG');
plot(ecg_filtered_isoline(1:n,1));
hold on; 
scatter(R_peaks_i(R_peaks_i<=n), ecg_filtered_isoline(R_peaks_i(R_peaks_i<=n),1), 'g', 'filled');
legend({'ECG signal', 'R peaks'});
xlabel('samples'); ylabel('voltage');
axis tight

% compute vectors that have two fiducial points for each complex/wave
Pwave_samples2= reshape(FPT_LeadI(:,[1,3])',1,[]);
QRS_samples2 = reshape(FPT_LeadI(:, [4,8])', 1, []);

figure; 
subplot(2,1,1)
title('Filtered ECG');
plot(ecg_filtered_isoline(1:n,1));
hold on; 
scatter(Pwave_samples2(Pwave_samples2<=n), ecg_filtered_isoline(Pwave_samples2(Pwave_samples2<=n),1), 'g', 'filled');
legend({'ECG signal','P wave'});

subplot(2,1,2)
plot(ecg_filtered_isoline(1:n,1));
hold on
scatter(QRS_samples2(QRS_samples2<=n), ecg_filtered_isoline(QRS_samples2(QRS_samples2<=n),1), 'r', 'filled');
xlabel('Samples'); ylabel('Amplitude (mV)');
legend({'ECG signal','QRS complex'});

%%% Compute segments time ranges
%RR distance
for i=2:length(R_peaks_i)
    dist_RR(i-1)=R_peaks_i(i)-R_peaks_i(i-1);
end
rangeRR=dist_RR/fs;
pt=10;
meanr_RR=trimmean(rangeRR,pt); %[s]

% %P wave duration
% for i=2:length(Pwave_samples2)
%     distP(i-1)= Pwave_samples2(i)-Pwave_samples2(i-1);
%   end
% distP(2:2:end)=[];
% rangeP=distP/fs;
% meanr_P=trimmean(rangeP,pt); %[s]

%QRS complex duration
for i=2:length(QRS_samples2)
    distQRS(i-1)=QRS_samples2(i)-QRS_samples2(i-1);
end
distQRS(2:2:end)=[];
rangeQRS=distQRS/fs;
meanr_QRS=trimmean(rangeQRS,pt); %[s]

%Identify the average duration of the PR segment, considering the distance between the first sample of the QRS complex and the first sample of the P wave
P_odd=Pwave_samples2(1:2:end);
QRS_odd=QRS_samples2(1:2:end);

for i=1:length(QRS_odd)
    PR_i(i)=QRS_odd(i)-P_odd(i);
end
rangePR=PR_i/fs;
meanr_PR=trimmean(rangePR,10);

%%% Trigger
% Initialize a variable to store msg
msg = {};

% WPW syndrome
if (round(meanr_PR,2) > 0.08 && round(meanr_PR,2) < 0.11) && round(meanr_QRS,2) >= 0.12 && round(meanr_RR,2) < 0.60
    msg{end+1} = 'Risk of WPW Syndrome.';
end

% first-degree AV block
if round(meanr_PR,2) >= 0.22 && round(meanr_QRS,2) < 0.12
    msg{end+1} = 'Risk of first-degree AV block.';
end

% first-degree AV block with wide QRS complex
if round(meanr_PR,2) >= 0.22 && round(meanr_QRS,2) > 0.12
    msg{end+1} = 'Risk of first-degree AV block with wide QRS complex.';
end

% second-degree AV block Mobitz type I 

% incomplete RBBB
if (round(meanr_QRS,2) >= 0.11 && round(meanr_QRS,2) < 0.12) 
    msg{end+1} = 'Risk of incomplete RBBB.';

% incomplete LBBB
elseif (round(meanr_QRS,2) >= 0.10 && round(meanr_QRS,2) < 0.12) 
    msg{end+1} = 'Risk of incomplete LBBB.';
end

% complete LBBB/RBBB
if round(meanr_QRS,2) >= 0.12 
    msg{end+1} = 'Risk of LBBB or RBBB.';
end

% tachycardia/bradycardia
if round(meanr_RR,2) < 0.60
    msg{end+1} = 'Abnormal RR interval: risk of tachycardia.';
end

if round(meanr_RR,2) > 1
    msg{end+1} = 'Abnormal RR interval: risk of bradycardia.';
end

% If no conditions were met
if isempty(msg)
    msg{end+1} = 'None of the tested pathologies has been detected.';
end

% Display all the msg
for i = 1:length(msg)
    warning(msg{i});
end

end





