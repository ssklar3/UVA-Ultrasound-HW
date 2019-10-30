ilam=2e-3;                       % wavelength
xs=[-7.5:.1:0]*1e-3;             % x field
zs=[0:.25:(25/(ilam*1e3))]*1e-3; % z field
y=0;                             % y value
frq=1e10;   % frequency of sampling
a=5e-3;     % x axis apature half width
b=7.5e-3;   % y axis apature half width
c=1540;     % speed of sound
% maxium time we want to calculate to
tmax=1.2*(sqrt((z^2)+((2*a+max(abs(x))).^2)+((2*b+max(abs(y))).^2))./c);
ps=zeros(length(xs),length(zs)); % pressure matrix initialization
%% function call
for n=1:length(zs)
    n
    for i=1:length(xs)
        [h,t,ps(i,n)]=rectapp(xs(i),y,zs(n),ilam,tmax);
    end
end
ps=ps./max(max(ps));    % normalize
%% plot
figure;
mesh((zs.*1e3)./(25./(ilam.*1e3)),abs((xs.*1e3)./5),ps)
xlabel('Depth (a^2/\lambda)')
ylabel('Azimuth (a)')
zlabel('Normalized Pressure')
title('Rectangular Mesh for a/\lambda=2.55')