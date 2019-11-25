%% Field call
path(path,'/Users/samuelsklar/Documents/Classwork/Ultrasound/Field_II_ver_3_24_mac');
field_init(-1);
% Set initial parameters
f0=3e6;               %  Transducer center frequency [Hz]
fs=100e6;             %  Sampling frequency [Hz]
c=1540;               %  Speed of sound [m/s]
lambda=c/f0;          %  Wavelength [m]
height=5/1000;
width=1/1000;
kerf=width/5;
N_elements=62;
N_active=32;
N_elements2=64;
focus=[0 0 50]/1000;  %  Initial electronic focus
N=200;
x_size = 20/1000;
y_size = 10/1000;
z_size = 20/1000;
z_start = 5/1000;
%  Number of scatterers
%  Width of phantom [mm]
%  Transverse width of phantom [mm]
%  Height of phantom [mm]
%  Start of phantom surface [mm];
%  Height of element [m]
%  Width of element [m]
%  Distance between transducer elements [m]
%  Number of elements
%  Number of elements
%  Define the transducers
Th = xdc_linear_array (N_elements, width, height, kerf, 1, 1, focus);
Th2 = xdc_linear_array (N_elements2, width, height, kerf, 1, 1, focus);
%  Set the impulse response and excitation of the emit aperture
impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (Th, impulse_response);
xdc_impulse (Th2, impulse_response);
excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (Th, excitation);
%  Define a small phantom with scatterers
%  Create the general scatterers
x = (rand (N,1)-0.5)*x_size;
y = (rand (N,1)-0.5)*y_size;
z = rand (N,1)*z_size + z_start;
positions=[x y z];
%  Generate the amplitudes with a Gaussian distribution
amp=randn(N,1);






[phantom_positions, phantom_amplitudes] = cyst_phantom(10000);
no_lines=N_elements-N_active+1; % Number of A-lines in image
dx=width; % Increment for image
z_focus=50/1000;
% Pre-allocate some storage
image_data=zeros(1,no_lines);
for i=1:no_lines
   
%     % Find position for imaging
%     x=(i-1-no_lines/2)*dx;
%     % Set the focus for this direction
%     xdc_center_focus (Th, [x 0 0]);
%     xdc_focus (Th, 0, [x 0 z_focus]);
%     xdc_center_focus (Th2, [x 0 0]);
%     xdc_focus (Th2, 0, [x 0 z_focus]);
%     % Set the active elements using the apodization
%     apo=[zeros(1, i-1) hamming(N_active)' zeros(1, N_elements-N_active-i+1)];
%     xdc_apodization (Th, 0, apo);
%     xdc_apodization (Th2, 0, apo);
    % Calculate the received response
    
%     [v, t1]=calc_scat(emit_aperture, receive_aperture, phantom_positions, phantom_amplitudes);
    
    [v, t1]=calc_scat_all(Th, Th2, phantom_positions, phantom_amplitudes,1);
    % Store the result
    image_data(1:(size(v,1)),1:(size(v,2)),i)=v;
    times(i) = t1;
end
% Free space for apertures
xdc_free (Th)
xdc_free (Th2)
% Adjust the data in time and display it as
% a gray scale image
min_sample=min(times)*fs;
for i=1:no_lines
    rf_env=abs(hilbert([zeros(round(times(i)*fs-min_sample),1); image_data(:,i)]));
    env(1:size(rf_env,1),i)=rf_env;
end
% make logarithmic compression to a 60 dB dynamic range
% with proper units on the axis
env_dB=20*log10(env);
env_dB=env_dB-max(max(env_dB));
env_gray=127*(env_dB+60)/60;
depth=((0:size(env,1)-1)+min_sample)/fs*c/2;
x=((1:no_lines)-no_lines/2)*dx;
image(x*1000, depth*1000, env_gray)
xlabel('Lateral distance [mm]')
ylabel('Depth [mm]')
axis('image')
colormap(gray(128))
title('Image of cyst phantom (60 dB dynamic range)')




%  Do the calculation
%[v,t]=calc_scat_all (Th, Th2, positions, amp, 1);
%  Plot the individual responses

% [N,M]=size(v);
% scale=max(max(v));
% v=v/scale;
% for i=1:M
%   plot((0:N-1)/fs+t,v(:,i)+i,'b'), hold on
% end
% hold off
% 
% title('Individual traces')
% xlabel('Time [s]')
% ylabel('Normalized response')
% axis([t t+N/fs 0 M+1])

