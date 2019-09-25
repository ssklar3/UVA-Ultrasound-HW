% Program to examine frequency dependent attenuation
% BME 4873/8730 use only
% HW Question: Extract the Gaussian function generation in the frequency domain.
% Simulate impact of tissue dependent attenuation  0.7 dB/cm/MHz
% Simulate and show the original spectrum (amplitude vs. frequency)
% for 80% -6dB fractional bandwidth (i.e. bandwidth between lower 
%     and upper frequency cutoffs at -6dB with respect to maximum 
%     divided by center frequency or mean of upper and lower cutoffs)
%     and 20% -6dB fractional bandwidth. Perform the calculation for
%     a 10 MHz center frequency and assume 3 cm depth of imaging.  
%     Remember to assume a two way path. 
%     Superimpose on the plot of the original spectrum the spectrum 
%     of the attenuated waveform and measure the center frequency 
%     downshift.

% (c) J A Hossack 2017,2019

clear all
%close all

% figure
fc=10e6;                                    % Center frequency
fs=fc/64;                                   % Define a sampling frequency (not critical)           
f=[fs:fs:2*fc];                             % Define an adequate frequency range
w=2*pi*f;                                   % Angular frequency radians
ns=length(f);                               % Number of samples

for i=1:2                                   % Case 20% then 80% BW
if i==1
    bw=20;                                  % Fractional bandwidth as percent
else
    bw=80;
end

sig=bw*fc/100;                              % Width of Gausian

gauss_pulse=exp(-pi*((f-fc)/sig).^2);       % Generate Gaussian pulse (frequency domain)

gauss_pulse=gauss_pulse./max(gauss_pulse);  % Normalize to 1 - although in this case it ought to be already
gauss_pulse_db=20*log10(gauss_pulse);

atten_fn_db=??????              % Atten using 0.7 db/MHz/cm assuming 2 x 3 cm pathlength (as an example)
                                % dB per MHz per cm. dB is "-" because it
                                % is attenuating
                                % f is f./1e-6 to account for f axis being
                                % in Hz and attenuation defined in terms of
                                % MHz
gauss_pulse_atten_db=gauss_pulse_db+atten_fn_db;  % Multiple in FD = add in log domain
figure
plot(f,gauss_pulse_db)
hold on
plot(f,(gauss_pulse_atten_db),'-');         % Attenuation function is a filter - multiplication in Freq. Dom. 
                                            % However, both filter and waveform are in log. domain so add them
plot(f,atten_fn_db,'--');
xlim([0 2e7]);
ylim([-80 0]);                              % Need to set useful limits
title('Original and attenuated spectrum')
xlabel('Frequency (Hz)')
ylabel('dB')
[Ymax,index] = max(gauss_pulse_atten_db);
f_max=f(index);
disp(['Center frequency of filtered function is ',num2str(f_max), ' Hz'])
% Hence calculate the frequency downshift from the original center freq.
% You can use data cursor in plot function to extract maximum if preferred
end