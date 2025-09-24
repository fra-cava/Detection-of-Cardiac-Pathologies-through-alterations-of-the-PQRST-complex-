clear all;
close all;
clc;

%WFDB Toolbox to load signal from Physionet 
wfdb2mat('100')
load('100m.mat');
ECG=val;
ECG=ECG.';
ecg = ECG;
[data, fs] = rdsamp('mitdb/100'); % specify sample rate

%the function needs two inputs: ecg and fs
[meanr_PR,meanr_QRS,meanr_RR,msg] = time_measures(ecg,fs);


%%
%load the signal that is acquired from BITalino
load('ECG_carly.mat');
ECG=signal;
ECG=ECG.';
ecg = ECG;
fs=1000;

%the function needs two inputs: ecg and fs
[meanr_PR,meanr_QRS,meanr_RR,msg] = time_measures(ecg,fs);
