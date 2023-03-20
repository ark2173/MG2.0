%% %Adapted from Robin Ji 2020
close all;
clear all;


File=['02-Jul-2022_Aged_WFM1_Sonication1'] %The name of the file without ]
addpath('/media/alina/Backup1/Programming/Matlab/Cavitation/')

%named
control = fread(fopen([File '_Control.bin'],'rb'),[338368, 100],'float32');
real = fread(fopen([File '_Real.bin'],'rb'),[338368, 1200],'float32');

% PCD analysis setup
PCD.harm = [3 4 5 6 7 8 9];% choose N-harmonic(s) to use
PCD.ultraharm = 3.5:8.5;
PCD.harmbw = 0.2;% bandwidth around harmonic
PCD.cavdose = [4.4 13.7].*1e6;% frequency range for cavitation dose
PCD.bb = 6.25e6;
PCD.bbbw = 0.1;

% time domain
Parameters.segmentsize = 505024;
Parameters.Fs= 50e6;
Parameters.freq = 1.5;
Parameters.dt = 1/Parameters.Fs;
Parameters.t = 0:Parameters.dt:(Parameters.segmentsize*Parameters.dt)-Parameters.dt;
% frequency domain
Parameters.df = (Parameters.freq*1e6)/(3000*6);
Parameters.Nfft = Parameters.Fs/Parameters.df;
Parameters.df = Parameters.Fs/Parameters.Nfft;
Parameters.f = 0:Parameters.df:Parameters.Fs-Parameters.df;
%Baseline data
[~,num_pulses]=size(control);
baseline = Filtering(control,PCD,num_pulses,Parameters);
%Signal data
[~,num_pulses]=size(real);
signals = Filtering(real,PCD,num_pulses,Parameters);

ICD = 10*log(rms(signals.bb-mean(baseline.bb))/rms(baseline.bb-mean(baseline.bb)))
SCD_UH = 10*log(rms(signals.ultraharm-mean(baseline.ultraharm))/rms(baseline.ultraharm-mean(baseline.ultraharm)))
SCD_H = 10*log(rms(signals.harm-mean(baseline.harm))/rms(baseline.harm-mean(baseline.harm)))

%Values
ICD=10*log(rms(signals.bb/rms(baseline.bb)));
SCD_UH=10*log(rms(signals.ultraharm/rms(baseline.ultraharm)));
SCD_H=10*log(rms(signals.harm/rms(baseline.harm)));