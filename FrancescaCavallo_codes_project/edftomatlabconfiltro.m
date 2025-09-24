
%script to convert a edf file (obtained from OpenSignals and Bitalino) to a
%variable in matlab , it must be adapted to your own file and variables.

clear; close all;
tt=edfread('ECG_Fra.edf') %tt is a kind of struct called "timetable", see its content; use the name of your edf file
fs=1000;%sample frequency, ensure from Opensignal
%signalA2=tt.SignalLabel2_RAW;%choose de channel Ax from tt (A2 in this case, which was used in bitalino to adquire, you have to adapt this line to your adquisition channel
signalA2=tt.RAW
nw=max(size(signalA2)); %signalA1 is formed by nw cells (including 1000 samples at each cell)
for i=1:nw 
    x(i,:)=signalA2{i,:}; %reading the cells and saving in a x matrix (nw x fm)
end
signal=reshape(x',1,nw*1000);%one raw whit the whole signal from x
%signal is ready for matlab .
%save('ECG_marta','fm','signal'); %to save the data in a .mat file, choose
%save('ECG_marta.mat','signal');
path='C:\Users\fraca\OneDrive\Desktop\Projects in biomedical engineering II\pan_tompkin\codes for ECG detection QRS\ECG_Processing'
save(fullfile(path, 'ECG_fra.mat'), 'signal');
%name and variables
%load('ECG_xxx.mat') ; %para load data from and existing mat file

%%
%example of representing in time
signaln=(signal-mean(signal))./std(signal); %normalized and center signal
nsamples=length(signaln);
Ts=1/fs; 
t=[0:nsamples-1]*Ts;
plot(t,-signaln);

%Filtering example, to eliminate the low frequencys below 1hz aprox (breathyng noise).
%om_a = 0.5 * 2*pi / fs;%discret=2pif/fs
%om_p = 2 * 2*pi / fs;
%Ra = 40;
%Rp = 1;
%transf bilineal
Wp = (2 / Ts) * tan(om_p / 2);
Wa = (2 / Ts) * tan(om_a / 2);
[orden, Wn] = buttord(Wp, Wa, Rp, Ra, 's');
[bs, as] = butter(orden, Wn, 'high', 's');

% Conversión del filtro prototipo del dominio continuo al discreto
[bz, az] = bilinear(bs, as, fs);
y=filter(bz,az,signaln);
figure(1), hold, plot(t, -y, 'r')