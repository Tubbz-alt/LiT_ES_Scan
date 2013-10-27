function plot_ps(beam,Nbin,fig_num,wait)

[nz,zb] = hist(beam(:,1),Nbin);
[ne,eb] = hist(beam(:,2),Nbin);
ps = hist2(beam(:,1),beam(:,2),zb,eb);

figure(fig_num);
subplot(2,2,1);
plot(ne,eb); axis tight;
subplot(2,2,2);
pcolor(zb,eb,ps'); shading flat;
subplot(2,2,4);
plot(zb,nz); axis tight;
if wait; pause; end