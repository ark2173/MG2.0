function baseline = Filtering(signals_control,PCD,NPulses,Parameters)
% Harmonic filter
harmfilt = zeros(Parameters.Nfft,1);
index_harm = (((Parameters.freq*1e6)/Parameters.df).*PCD.harm)+1;
harm_bw = (PCD.harmbw*1e6)/(2*Parameters.df);
for ind = 1:length(PCD.harm)
    harmfilt(index_harm(ind)-harm_bw:index_harm(ind)+harm_bw) = 1;
end

% UltraHarmonic filter
ultraharmfilt = zeros(Parameters.Nfft,1);
index_ultraharm = (((Parameters.freq*1e6)/Parameters.df).*PCD.ultraharm)+1;
harm_bw = (PCD.harmbw*1e6)/(2*Parameters.df);
for ind = 1:length(PCD.ultraharm)
    ultraharmfilt(index_ultraharm(ind)-harm_bw:index_ultraharm(ind)+harm_bw) = 1;
end

% Cavitation dose filter
dosefilt = zeros(Parameters.Nfft,1);
index_dose = PCD.cavdose/Parameters.df;
dosefilt(index_dose(1):index_dose(2)) = 1;

% Broadband filter
bbfilt = zeros(Parameters.Nfft,1);
index_bb = PCD.bb/Parameters.df;
bb_bw = (PCD.bbbw*1e6)/(2*Parameters.df);
bbfilt(index_bb-bb_bw:index_bb+bb_bw) = 1; 

control_dosefilt = repmat(dosefilt,1,NPulses);
control_harmfilt = repmat(harmfilt,1,NPulses);
control_ultraharmfilt = repmat(ultraharmfilt,1,NPulses);
control_bbfilt = repmat(bbfilt,1,NPulses);
curr_control_sig_fft = fft(signals_control,Parameters.Nfft,1);

baseline.dose = sum(abs(curr_control_sig_fft.*control_dosefilt).^2);
baseline.harm = sum(abs(curr_control_sig_fft.*control_harmfilt).^2);
baseline.ultraharm = sum(abs(curr_control_sig_fft.*control_ultraharmfilt).^2);
baseline.bb = sum(abs(curr_control_sig_fft.*control_bbfilt).^2);
baseline.mean_harm = mean(baseline.harm);
baseline.mean_ultraharm = mean(baseline.ultraharm);
baseline.mean_bb = mean(baseline.bb);
               
end