%% Field call
path(path,'/Users/samuelsklar/Documents/Classwork/Ultrasound/Field_II_ver_3_24_mac');
field_init(-1);
%% Variables
a=5e-3;                 % Half x width
b=7.5e-3;               % Half Y width
lam=2e-3;               % Wavelength
z=(.1:.1:12.5).*1e-3;   % Z field
y=(0).*7.5e-3;          % Y field
x=(0:.1:1.5).*5e-3;     % X field
xpitch=.04e-3;          % Element sizes in x dimension
ypitch=.04e-3;          % Element sizes in y dimension
xpts=-a:xpitch:a;        % X co ordinates of element corners
ypts=-b:ypitch:b;        % Y coordinates of element corners
noels=(length(xpts)-1)*(length(ypts)-1); % Number of elements
%% Rectangles input definition
rect=zeros(noels,19);
rect(:,1)=1;
rect(:,14)=1;
rect(:,15)=xpitch;
rect(:,16)=ypitch;
count=0;
for j=1:length(ypts)-1
    for k=1:length(xpts)-1
        count=count+1;
        rect(count,2:3)=[xpts(k),ypts(j)];
        rect(count,5:6)=[xpts(k+1),ypts(j)];
        rect(count,8:9)=[xpts(k+1),ypts(j+1)];
        rect(count,11:12)=[xpts(k),ypts(j+1)];
        rect(count,17:18)=[mean(xpts(k:k+1)),mean(ypts(j:j+1))];
    end
end
%% Field Parameters
frq=1e10;
set_sampling(frq);
center=[0,0,0];
focus=[0,0,12.5e-3];
tmax=4e-5;
c=1540;
%% Field call and pressure calc
Th = xdc_rectangles (rect, center, focus);  % Transducer definition

p=zeros(length(x),length(z));            % Pressure array initialization
for i=1:length(z)
    i
    for n= 1:length(x)
        [h,t]=calc_h (Th,[x(n),y,z(i)]); % Transfer function calc
        t=[t:1/frq:tmax];                % Define times for sin wave
        wav=sin(2*pi*t*c./lam);          % creation of sin wave
        hwc=conv(h,wav);                 % convolution
        hwc=diff(hwc);                   % differentiation
        prS=floor((length(hwc)./2)-(length(hwc)./5));   % Lower bounding for amplitude check
        prE=ceil((length(hwc)./2)+(length(hwc)./5));    % Higher bounding for amplitude check
        p(n,i)=max(hwc(prS:prE))-min(hwc(prS:prE));          % Calculate pressure
    end
end
p=p./max(max(p)); % Pressure Normalization
%% plot
figure
    mesh(z/((a^2)/lam),x/(5e-3),p);
    xlabel('Depth (a^2/\lambda)')
    ylabel('Azimuth (a)')
    zlabel('Normalized Pressure')
    title('Rectangular Mesh for a/\lambda=2.5')