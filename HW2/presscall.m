ilam=1e-3;                       % wavelength
xs=[-7.5:.1:0]*1e-3;             % x field
zs=[0:.25:(25/(ilam*1e3))]*1e-3; % z field
ps=zeros(length(xs),length(zs)); % pressure matrix initialization
for n=1:length(zs)
    n
    for i=1:length(xs)
        ps(i,n)= pressfun(zs(n),xs(i),ilam); % call pressure calculation fn
    end
end
ps=ps./max(max(ps));    % normalize
%% plot
figure;
mesh((zs.*1e3)./(25./(ilam.*1e3)),abs((xs.*1e3)./5),ps)
xlabel('Depth (a^2/\lambda)')
ylabel('Azimuth (a)')
zlabel('Normalized Pressure')
title('Circular Mesh for a/\lambda=5')