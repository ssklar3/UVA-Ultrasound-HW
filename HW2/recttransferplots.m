x=[1.4,1,0.6,0.6,0].*5e-3;  % Set of x values
y=[1.25,0.75,0.75,0,0].*7.5e-3; % Set of y values
z=50e-3; % Z value
lam=1e-3; % a lam just so the code will run
frq=1e10; % frequency of sampling
a=5e-3;   % x dimension half width
b=7.5e-3; % y dimension half width
c=1540;   % speed of sound in body
% maxium t we want to calculate values to
tmax=1.2*(sqrt((z^2)+((2*a+max(abs(x))).^2)+((2*b+max(abs(y))).^2))./c);
t=zeros(length(x),floor(tmax*frq)); % bin for times
h=t; % bin for transfer function values
p=zeros(length(x),1); % bin for pressures that we dont use here
%% plot
figure;
hold on
for i=1:length(x)
    [h(i,:),t(i,:),p(i,:)]=rectapp(x(i),y(i),z,lam,tmax); % call function
    subplot(length(x),1,i)
    plot(t(i,:),h(i,:))
    xlabel('Time (s)')
    ylabel('h(x,t)/c')
    axis([3.2e-5 3.5e-5 0 1.2])
    title(['x/a = ',num2str(x(i)/5e-3),'y/b = ',num2str(y(i)/7.5e-3)])
end