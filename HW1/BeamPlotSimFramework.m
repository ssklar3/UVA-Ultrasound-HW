% Program to simulate beam response using Huygen's Principle
% BME 8730 use only
% (c) J A Hossack 2019
% MKS: meters kilograms seconds - recommended to use this to avoid losing
% track
later (using mm would typically be more convenient for this scale)

% For the questions, you change center frequency (fc), pitch, num_elems,
% focal depth and steer angle. You should be able to find and change these
% without difficulty
% You need to enter fresh code in the key missing lines using discussion
% from class and use of geometry to calculation focal and propagation
% delays

% Propagation delays - actual time of flight from individual element to
% point in field in question (i.e. we loop over a range of points in the
% field) Propagation delays change for every point in the field examined

% Focal delays compensatory delays applied to compensate for actual
% propagration delays. Focal delays are a function of array geometry and
% focal requirement. They are fixed once they are calculated at the
% beginning of the problem.
% In a real system, we are concerned with relative focal delays - i.e. the
% minimum focal delay (or time advance) is substracted from all calculated
% focal delays (so that the minimum time delay is zero and not some needlessly long delay
% We do not need this step in this simulation so it is
% omitted entirely. We discuss this in class. Dont worry about this detail.

% At the focus, all propagation delays and all focal delays cancel for all
% array elements. Since focal delays are compensatory, we think of them as
% having "negative" delays. i.e. a negative time delay (time advance) is
% used so that then added to the positive propagation delays, apparent
% delay precisely cancels.

clear all                                   % Clear work space from previous runs etc
%close all                                  % Close old figures

% figure

% Be sure to set all parameters for each problem according to the problem
z_foc=30e-3;                                % Range direction focal distance
theta_steer=0*pi/180;                       % Steer angle - will use this later


num_elems=128;                              % Number of array elements

   
fc=10e6;                                    % Center frequency
pitch=0.075e-3;                             % Array center-center element spacing
                                            % Set pitch to multiples of
                                            % this to test impact of
                                            % undersampling

%weight=hann(num_elems);                    % Apodization / weighting
                                            %function
weight=ones(num_elems,1);  

        
vel=1540;                                   % Speed of sound - all units MKS

fs=fc/64;                                   % Define a sampling frequency (not critical)
                                            % Observe result when we choose
                                            % "64" to be "32" or "128".
                                            % This is sampling in the
                                            % frequency domain so it does
                                            % what when we convert to time
                                            % domain?
                                            
f=[fs:fs:8*fc];                             % Define an adequate frequency range (improve resolution in time domain)
                                            % Observe the result when you
                                            % make the "8" a "4" or a "16"
                                            % -Plot spectrum and plot time
                                            % domain - essentially, the
                                            % maximum frequency is extended
                                            % and this is the Nyquist limit
                                            % so it does what in the time
                                            % domain?
i_cen=round(fc/fs);                          % Index of center frequency for single frequency 
                                            % calculation

w=2*pi*f;                                   % Angular frequency radians
ns=length(f);                               % Number of samples

tdel=0e-6;                                 % Use a fixed time offset so that initial waveform is 0 for t<0 (not critical)
tdel2=1.0e-6;                               % Use a fixed time offset so that initial waveform is 0 for t<0 (not critical)
                                            % Plot the time domain version
                                            % using each of these delays
                                            % and observe the effect

bw=30;                                      % Fractional bandwidth as percent
bw2=80;                                     % Fractional bandwidth as percent
sig=bw*fc/100;                              % Width of Gausian
gauss_pulse=exp(-pi*((f-fc)/sig).^2);       % Generate Gaussian pulse (frequency domain)

gauss_pulse=gauss_pulse.*exp(-j*w*tdel);    % Apply time delay so 0 for t<0

figure
plot(f,20*log10(gauss_pulse));
ylim([-60,0])
% xlim([0,2e7])                                 
title('Gaussian input function spectrum')
xlabel('Frequency (Hz)')
ylabel('dB')
f_max=f(ns)                                 % Just a verification readout of maximum frequency


gauss_t=real(ifft(gauss_pulse));            % Time domain of base waveform for reference
gauss_t=gauss_t./max(gauss_t);
env_gauss_t=abs(hilbert(gauss_t));          % Envelope calculation - i.e. peak value of sinusoidal function

tstep=1./max(f);                            % Time steps after using Inverst FFT
t=[1:ns].*tstep;                            % Define time axis
% 
figure                                    % Plot waveform as a diagnostic
plot(t,gauss_t)
hold on
plot(t,env_gauss_t,'--')
title('Gaussian Pulse')
xlabel('Time (s)')
ylabel('Pressure - Arb Units')

%weight=ones(num_elems,1);                  % Define the weighting function (you can change this)
%weight=hann(num_elems);


% Be sure to reset the sampling according to the problem so that you
% achieve a smooth reconstruction of the main beam plus at least a few
% sidelobes if present. You may want to run a few times and determine best
% values. To the extent you can use similar sampling to make comparisons
% easier between successive questions, do so. It will help answer
% questions.
% You can  read out from figures using Matlab data cursor functions under
% tools or cross icon.

if 0
x_pts=[-5:0.2:5].*1e-3;                     % Define X-direction field locations
z_pts=[5.01:1.0:50.01].*1e-3;               % Define Z-direction field locations

% Create focal delays
for i=1:num_elems
    x_elem(i)=((i-1)-(num_elems-1)./2).*pitch;  % Calculate locations of each array element

    
    %foc_del(i) = ??????
    %steer_del(i) = ?????? (when steering delay is applied)
                                         
end

for j=1:length(z_pts)                       % Loop over field locations
    fprintf('%d ',j)                        % Short piece of code to provide record of progress
    if (j/20)==round(j/20)
            fprintf('\n')
    end
    

    for i=1:length(x_pts)
        sum_pulse=zeros(size(gauss_pulse));                 % Initialize sum of waveforms to zero before starting loop
        for k=1:num_elems

            
            %prop_del = ???? for each elem location, calculate actual
           
            % propagation time out to current field point
            sum_pulse=sum_pulse+weight(k).*gauss_pulse.*... % Sum pulse taking account of current pulse, propagation delay - focusing delay
                exp(-sqrt(-1)*w.*(prop_del-foc_del(k)));
            
        end

        sum_pulse_t=real(ifft(sum_pulse));                  % Convert to time domain

        field_val(j,i)= max(abs(hilbert(sum_pulse_t)));     % Field values in beamplot - Hilbert finds envelope of sinsoidal fn
        field_val_fc(j,i)=abs(sum_pulse(i_cen));            % Magnitude of summed single frequency (center frequency)
 
    end
    field_val(j,:)=field_val(j,:)./max(field_val(j,:));     % Normalize by peak value at each range
    field_db(j,:)=20*log10(field_val(j,:));                 % Convert to dB
    field_val_fc(j,:)=field_val_fc(j,:)./max(field_val_fc(j,:));
    field_fc_db(j,:)=20*log10(field_val_fc(j,:)); 
end

% Note - in the above normalization, the normalization is performed once
% per range location. This is common in beamplots related to practical
% scanners. However, when calculating true 2D beamplot it is sometimes
% better to normalize by the peak value across the entire field and not on
% a per range basis.

figure
contour(z_pts*1e3,x_pts*1e3,field_db',[-30 -20 -15 -12 -9 -6 -3]) % Contour plot, defined contours
%mesh(z_pts*1e3,x_pts*1e3,field_db');
%plot(x_pts*1e3,field_db)
%axis([min(x_pts*1e3) max(x_pts*1e3) -40 0])
%colormap('gray')
%hold on
title(title_label)                          % Label according to question
xlabel('Range (mm)')                        % Be sure to account for mm or m here
ylabel('Azimuth (mm)')
axis equal                                  % Makes the X and Y indexing match (correct aspect ratio)

end