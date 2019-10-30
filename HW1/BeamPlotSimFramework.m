z_foc=50e-3;                                % Range direction focal distance\
theta_steer=(45)*pi/180;
num_elems=64;                              % Number of array elements
fc=10e6;                                    % Center frequency
pitch=(.075e-3);                         
vel=1540;                                   % Speed of sound - all units MKS
fs=fc/64;                                   % Define a sampling frequency (not critical)                                         
f=[fs:fs:8*fc];                             % Define an adequate frequency range (improve resolution in time domain)
i_cen=round(fc/fs);                         % Index of center frequency for single frequency calc
w=2*pi*f;                                   % Angular frequency radians
ns=length(f);                               % Number of samples
tdel=0e-6;                                  % Use a fixed time offset so that initial waveform is 0 for t<0 (not critical)
tdel2=1.0e-6;                               % Use a fixed time offset so that initial waveform is 0 for t<0 (not critical)
bw=30;                                      % Fractional bandwidth as percent
bw2=80;                                     % Fractional bandwidth as percent
sig=bw*fc/100;                              % Width of Gausian
gauss_pulse=exp(-pi*((f-fc)/sig).^2);       % Generate Gaussian pulse (frequency domain)
gauss_pulse=gauss_pulse.*exp(-j*w*tdel);    % Apply time delay so 0 for t<0
%att= .00005*f;
att=0;
%%
weight=ones(num_elems,1);                  % Define the weighting function (you can change this)

%weight=hann(num_elems);
% x_pts=[-2:0.2:2].*1e-3;    
% z_pts=[5:1.0:2*z_foc*10^3].*1e-3;
x_pts=[-60:0.2:60].*1e-3;                     % Define X-direction field locations
z_pts=[5:1.0:60].*1e-3;

x_elem=zeros(num_elems,1);
foc_del=zeros(num_elems,1);
steer_del=zeros(num_elems,1);
for i=1:num_elems
    x_elem(i)=((i-1)-(num_elems-1)./2).*pitch;  % Calculate locations of each array element]
    steer_del(i) =  (((pitch*sin(theta_steer)*i)+sqrt(((x_elem(i))^2)+(z_foc^2))))./vel;  
end
field_val= zeros(length(z_pts),length(x_pts));
field_val_fc= zeros(length(z_pts),length(x_pts));
field_db= zeros(length(z_pts),length(x_pts));
field_fc_db= zeros(length(z_pts),length(x_pts));
for j=1:length(z_pts)                       % Loop over field locations
    fprintf('%d ',j)                        % Short piece of code to provide record of progress
    if (j/20)==round(j/20)
            fprintf('\n')
    end
    for i=1:length(x_pts)
        sum_pulse=zeros(size(gauss_pulse));                 % Initialize sum of waveforms to zero before starting loop
        for k=1:num_elems
            r=sqrt((x_elem(k)-x_pts(i))^2+z_pts(j)^2);
            prop_del = (r./vel);
              sum_pulse=sum_pulse+weight(k).*10.^(-att.*r./20).*gauss_pulse.*... % Sum pulse taking account of current pulse, propagation delay - focusing delay
                exp(-sqrt(-1)*w.*(prop_del-steer_del(k)));
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
%%
figure;
contour(z_pts*1e3,x_pts*1e3,field_db',[-30 -20 -15 -12 -9 -6 -3]);
title('Contour Field')
xlabel('Range (mm)')
ylabel('Azimuth (mm)')
figure
mesh(z_pts*1e3,x_pts*1e3,field_db');
title('Mesh Field')
xlabel('Range (mm)')
ylabel('Azimuth (mm)')
zlabel('Signal (dB)')
figure
plot(x_pts*1e3,field_db(-5+z_foc*10^3,:))
title('Beam at Focus')
ylabel('Signal (dB)')
xlabel('Azimuth (mm)')