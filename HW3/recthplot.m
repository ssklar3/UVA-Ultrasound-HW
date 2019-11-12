%% Field call
path(path,'/Users/samuelsklar/Documents/Classwork/Ultrasound/Field_II_ver_3_24_mac');
field_init(-1);
%% Variables
a=5e-3;                 % Half x width
b=7.5e-3;               % Half Y width
y=[1.25,0.75,0.75,0,0].*7.5e-3;          % Y points
x=[1.4,1,0.6,0.6,0].*5e-3;               % X points
z=(50).*1e-3;           % Z point
z=z*ones(1,length(x));
pitch=.02e-3;           % Element sizes
xpts=-a:pitch:a;        % X coordinates of element corners
ypts=-b:pitch:b;        % Y coordinates of element corners
noels=(length(xpts)-1)*(length(ypts)-1); % Number of elements
%% Rectangles input definition
rect=zeros(noels,19);
rect(:,1)=1;
rect(:,14)=1;
rect(:,15)=pitch;
rect(:,16)=pitch;
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
frq=1e11;
set_sampling(frq);
center=[0,0,0];
focus=[0,0,50e-3];
tmax=4e-5;
c=1540;
%% Field call
Th = xdc_rectangles (rect, center, focus);  % Transducer definition
[h,t]=calc_h (Th,[x',y',z']); % Transfer function calc
h=h./max(max(h));
%% Add padding on front end
hfront=zeros(ceil(((t-(3e-5))*frq)),5);
tf=(3e-5:1/frq:t+length(h)/frq);
h=[hfront;h];
%% plot
figure
for i=1:length(x)
    subplot(length(x),1,i)
    plot(tf,h(:,i));
    xlabel('Time (s)')
    ylabel('h(x,t)/c')
    axis([3.2e-5 3.6e-5 0 1.2])
    title(['x/a = ',num2str(x(i)/5e-3),'y/b = ',num2str(y(i)/7.5e-3)])
end

    