
                                            % Observe result when we choose
                                            % "64" to be "32" or "128".
                                            % This is sampling in the
                                            % frequency domain so it does
                                            % what when we convert to time
                                            % domain?

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