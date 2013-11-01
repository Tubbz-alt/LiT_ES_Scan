function plot_ps(LiT_struct,Nbin,normalize,fig_num,wait,use_data,yag_axis,yag_spec,tcav_axis,tcav_spec)

beam = LiT_struct.BEAM;
if use_data == 1
    [beam, ~, ~, yag_axis, ~, yag_spec] = comp_data_litrack(beam,Nbin,yag_axis,yag_spec);
    [nz,zb] = hist(beam(:,1),Nbin);
    [ne,eb] = hist(beam(:,2),yag_axis);
    ne(1) = 0;
    ne(end) = 0;
    ps = hist2(beam(:,1),beam(:,2),zb,eb);
elseif use_data == 2
    [beam, tcav_axis, ~, yag_axis, ~, yag_spec, tcav_spec] = comp_data_litrack(beam,Nbin,yag_axis,yag_spec,tcav_axis,tcav_spec);
    [nz,zb] = hist(beam(:,1),tcav_axis);
    nz(1) = 0;
    nz(end) = 0;
    [ne,eb] = hist(beam(:,2),yag_axis);
    ne(1) = 0;
    ne(end) = 0;
    ps = hist2(beam(:,1),beam(:,2),zb,eb);
elseif normalize
    [z_proj,e_proj,ps] = gen_ps(LiT_struct,Nbin);
    zb = z_proj(:,1);
    nz = z_proj(:,2);
    eb = e_proj(:,1);
    ne = e_proj(:,2);
else
    [nz,zb] = hist(beam(:,1),Nbin);
    [ne,eb] = hist(beam(:,2),Nbin);
    ps = hist2(beam(:,1),beam(:,2),zb,eb);
end

figure(fig_num);

subplot(2,2,1);
plot(ne,eb,'b','linewidth',2);
if normalize || use_data 
    ylabel('\delta [%]','fontsize',14);
    if normalize
        xlabel('N/bin','fontsize',14);
    end
else
    ylabel('Energy [GeV]','fontsize',14);
end
if use_data
    hold on;
    plot(yag_spec,yag_axis,'r','linewidth',2);
    hold off;
end
axis tight;

subplot(2,2,2);
pcolor(zb,eb,ps'); shading flat;
if normalize || use_data 
    ylabel('\delta [%]','fontsize',14);
    xlabel('Z [\mum]','fontsize',14);
else
    ylabel('Energy [GeV]','fontsize',14);
    xlabel('Z [m]','fontsize',14);
end

subplot(2,2,4);
plot(zb,nz,'b','linewidth',2);
if normalize || use_data == 2
    xlabel('Z [\mum]','fontsize',14);
    if normalize
        ylabel('I [kA]','fontsize',14);
    end
else
    xlabel('Z [m]','fontsize',14);
end
if use_data == 2
    hold on;
    plot(tcav_axis,tcav_spec,'r','linewidth',2);
    hold off;
end
axis tight;

if wait; pause; end