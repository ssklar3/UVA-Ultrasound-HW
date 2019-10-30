bws=[30,50,80];
fcs=[2e6,10e6];
z_foc=[1e-2,3e-2];
GAFout=cell(length(bws),length(fcs),length(z_foc));
attenuation=zeros(length(bws),length(fcs),length(z_foc));
acf=zeros(length(bws),length(fcs),length(z_foc));
count=0;
figure
for m=1:length(bws)
    
    for n=1:length(fcs)
        
       for k=1:length(z_foc)
           count = count+1;
           [f,att,gauss_pulse,gauss_pulse_att,gauss_t,env_gauss_t,t]=GauAttFn(bws(m),fcs(n),z_foc(k),64,8);
           attenuation(m,n,k) = max(20*log10(gauss_pulse_att));
           acf(m,n,k) = spline(20*log10(gauss_pulse_att(20:64)),f(20:64),attenuation(m,n,k));
            subplot(3,4,count);
            plot(f,20*log10(gauss_pulse));
            hold on
            plot(f,att);
            plot(f,20*log10(gauss_pulse_att));
            plot(f,20*log10(gauss_pulse_att)-attenuation(m,n,k),'--');
            ylim([-60,0])
            xlim([fcs(n)*0,fcs(n)*2])                                 
            title([num2str(fcs(n)*10^-6),' MHz, ', num2str(z_foc(k)*100),...
                ' cm, -6dB BW ', num2str(bws(m)), '%, Att:',...
                num2str(round(attenuation(m,n,k),2)), 'dB, ACF: ',...
                num2str(round(acf(m,n,k)*10^-6,2)),'MHz'])
            xlabel('Frequency (Hz)')
            ylabel('dB')
            if count==1
                legend('Original Gaussian','Attenuation','Attenuated Gaussian'...
                    ,'Renormalized Attenuated Gaussian')
            end
            hold off
       end
       
    end
end