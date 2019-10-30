z_foc=50e-3;                                % Range direction focal distance
theta_steer=(0)*pi/180;
num_elems=32;                              % Number of array elements
fc=5e6;                                    % Center frequency
pitch=(.075e-3)*2;                         
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
%% set x_pts to [0] for axis view

x_pts=[-10:0.2:10].*1e-3;
%y_pts=[-5:0.2:5].*1e-3;
y_pts=0;
%x_pts=0;
z_pts=[2:.2:80].*1e-3;
x_elem=[-5:0.5:5].*1e-3;
y_elem=[-5:0.5:5].*1e-3;
w_elem=zeros(length(x_elem),length(y_elem));
foc_del=zeros(num_elems);
for i=1:length(x_elem)
    %x_elem(i)=((i-1)-(num_elems-1)./2).*pitch;  % Calculate locations of each array element]
    for k=1:length(y_elem)
        %y_elem(k)=((k-1)-(num_elems-1)./2).*pitch;
        foc_del(i,k) =  (sqrt(((x_elem(i))^2)+(y_elem(k)^2)+(z_foc^2)))./vel;
        if sqrt(((x_elem(i))^2)+(y_elem(k)^2))<=5e-3
            w_elem(i,k)=1;
        end
    end
end

weight=ones(length(x_elem),length(y_elem));
%weight=hann(length(x_elem))*hann(length(y_elem))';

field_val= zeros(length(z_pts),length(x_pts),length(y_pts));
field_val_fc= zeros(length(z_pts),length(x_pts),length(y_pts));
field_db= zeros(length(z_pts),length(x_pts),length(y_pts));
field_fc_db= zeros(length(z_pts),length(x_pts),length(y_pts));
for q=1:length(z_pts)
        fprintf('%d ',q)                        % Short piece of code to provide record of progress
    if (q/20)==round(q/20)
            fprintf('\n')
    end
    for i=1:length(x_pts)
        for n=1:length(y_pts)
        sum_pulse=zeros(size(gauss_pulse));                 % Initialize sum of waveforms to zero before starting loop
            for k=1:length(y_elem)
                for m=1:length(y_elem)
                r=sqrt(((x_elem(k)-x_pts(i))^2)+((y_elem(m)-y_pts(n))^2)+(z_pts(q)^2));
                prop_del = (r./vel);
                  sum_pulse=sum_pulse+(w_elem(k,m)*weight(k,m)*r^-2).*10.^(-att.*r./20).*gauss_pulse.*... % Sum pulse taking account of current pulse, propagation delay - focusing delay
                    exp(-sqrt(-1)*w.*(prop_del-foc_del(k,m)));
%                   sum_pulse=sum_pulse+(w_elem(k,m)*weight(k,m)).*10.^(-att.*r./20).*gauss_pulse.*... % Sum pulse taking account of current pulse, propagation delay - focusing delay
%                     exp(-sqrt(-1)*w.*(prop_del-foc_del(k,m)));
                end
            end
        
        sum_pulse_t=real(ifft(sum_pulse));                  % Convert to time domain

        field_val(q,i,n)= max(abs(hilbert(sum_pulse_t)));     % Field values in beamplot - Hilbert finds envelope of sinsoidal fn
        field_val_fc(q,i,n)=abs(sum_pulse(i_cen));            % Magnitude of summed single frequency (center frequency)
        end
    end
% Typical normalization
      field_val(q,:,:)=field_val(q,:,:)./max(max(field_val(q,:,:)));     % Normalize by peak value at each range
      field_db(q,:,:)=20*log10(field_val(q,:,:));                 % Convert to dB
      field_val_fc(q,:,:)=field_val_fc(q,:,:)./max(max(field_val_fc(q,:,:)));
      field_fc_db(q,:,:)=20*log10(field_val_fc(q,:,:)); 
end
%% Axis FC normalization
%       field_val=field_val./max(max(max(field_val)));     % Normalize by peak value at each range
%       field_db(q,:,:)=20*log10(field_val(q,:,:));                 % Convert to dB
%       field_val_fc(:,:)=field_val_fc(:,:)./max(max(field_val_fc(:,:)));
%       field_fc_db=20*log10(field_val_fc); 
%%
figure;
contour(x_pts*1e3,z_pts*1e3,squeeze(field_db(:,:)),[-30 -20 -15 -12 -9 -6 -3]);
%contour(squeeze(field_db(:,:,26)),[-30 -20 -15 -12 -9 -6 -3]);
title('Contour Field')
ylabel('Range (mm)')
xlabel('Azimuth (mm)')
figure
mesh(x_pts*1e3,z_pts*1e3,squeeze(field_db(:,:)));
%mesh(squeeze(field_db(:,:,26)));
title('Mesh Field')
ylabel('Range (mm)')
xlabel('Azimuth (mm)')
zlabel('Signal (dB)')
figure
mesh(x_pts*1e3,z_pts*1e3,squeeze(field_fc_db(:,:)));
%mesh(squeeze(field_db(:,:,26)));
title('Mesh Field at FC')
ylabel('Range (mm)')
xlabel('Azimuth (mm)')
zlabel('Signal (dB)')
% figure
% plot(z_pts*1e3,squeeze(field_fc_db(:,:)));
% title('Center axis at FC')
% xlabel('Range (mm)')
% ylabel('Signal (dB)')