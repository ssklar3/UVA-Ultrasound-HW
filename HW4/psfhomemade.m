path(path,'/Users/samuelsklar/Documents/Classwork/Ultrasound/Field_II_ver_3_24_mac');
field_init(-1);
%% 
bw=30;                                      % bandwidth
fc=10e6;
z_foc= 50e-3;
sf=[64];
range=[8];
fs=fc/sf;                                   % Define a sampling frequency (not critical)
f=[fs:fs:range*fc];                         % Define an adequate frequency range (improve resolution in time domain)
w=2*pi*f;                                   % Angular frequency radians
ns=length(f);                               % Number of samples
tdel=0e-6;                                  % Use a fixed time offset so that initial waveform is 0 for t<0 (not critical)
tdel2=1.0e-6;                               % Use a fixed time offset so that initial waveform is 0 for t<0 (not critical)
                                            % Plot the time domain version
                                            % using each of these delays
                                            % and observe the effect
sig=bw*fc/100;                              % Width of Gausian
gauss_pulse=exp(-pi*((f-fc)/sig).^2);       % Generate Gaussian pulse (frequency domain)

gauss_pulse=gauss_pulse.*exp(-j*w*tdel);    % Apply time delay so 0 for t<0
att= -.00005*z_foc*f*2;

gauss_pulse_att=10.^(att./20).*gauss_pulse;         % Attenuation
psf=ifft(gauss_pulse);            % Time domain of base waveform for reference
%%

f0=3e6;                  %  Transducer center frequency [Hz]
fs=100e6;                %  Sampling frequency [Hz]
c=1540;                  %  Speed of sound [m/s]
lambda=c/f0;             %  Wave length [m]
width=lambda;            %  Width of element
element_height=5/1000;   %  Height of element [m]
kerf=width/20;           %  Kerf [m]
focus=[0 0 50]/1000;     %  Fixed focal point [m]
N_elements=192;          %  Number of elements in the transducer
N_active=64;             %  Active elements in the transducer
%  Set the sampling frequency
set_sampling(fs);
%  Generate aperture for emission
emit_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 5, focus);
%  Set the impulse response and excitation of the emit aperture
impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (emit_aperture, impulse_response);
excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (emit_aperture, excitation);
%  Generate aperture for reception
receive_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 5, focus);
%  Set the impulse response for the receive aperture
xdc_impulse (receive_aperture, impulse_response);
%   Load the computer phantom
[phantom_positions, phantom_amplitudes] = cyst_phantom(10000);