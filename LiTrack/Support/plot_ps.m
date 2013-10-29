function plot_ps(beam,Nbin,fig_num,wait,yag_axis,yag_spec,tcav_axis,tcav_spec)

lit_z = 1e6*beam(:,1);
lit_d = 100*(beam(:,2)-mean(beam(:,2)))/mean(beam(:,2));

[nz,zb] = hist(lit_z,Nbin);
[ne,eb] = hist(lit_d,Nbin);

if nargin == 8
    
    [~,lit_ind] = max(nz);
    lit_z = lit_z - zb(lit_ind);

    [~,dat_ind] = max(tcav_spec);
    tcav_axis = tcav_axis - tcav_axis(dat_ind);
    
    [nz,zb] = hist(lit_z,tcav_axis);
    nz(1) = 0;
    nz(end) = 0;
    lit_zmax = max(nz);
    dat_zmax = max(tcav_spec);
    nz = (dat_zmax/lit_zmax)*nz;

    [ne,eb] = hist(lit_d,yag_axis);
    ne(1) = 0;
    ne(end) = 0;
    lit_emax = max(ne);
    dat_emax = max(yag_spec);
    ne = (dat_emax/lit_emax)*ne;

elseif nargin == 6

    [nz,zb] = hist(lit_z,Nbin);
    [ne,eb] = hist(lit_d,yag_axis);
    ne(1) = 0;
    ne(end) = 0;
    lit_emax = max(ne);
    dat_emax = max(yag_spec);
    ne = (dat_emax/lit_emax)*ne;

end

ps = hist2(lit_z,lit_d,zb,eb);

figure(fig_num);

subplot(2,2,1);
plot(ne,eb,'b','linewidth',2); 
if nargin >= 6
    hold on;
    plot(yag_spec,yag_axis,'r','linewidth',2);
    hold off;
end
axis tight;
subplot(2,2,2);
pcolor(zb,eb,ps'); shading flat;
subplot(2,2,4);
plot(zb,nz,'b','linewidth',2);
if nargin >= 8
    hold on;
    plot(tcav_axis,tcav_spec,'r','linewidth',2);
    hold off;
end
axis tight;

if wait; pause; end