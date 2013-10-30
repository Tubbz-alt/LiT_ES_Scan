function plot_ps(beam,Nbin,fig_num,wait,yag_axis,yag_spec,tcav_axis,tcav_spec)

if nargin == 8
    [nz,zb] = hist(beam(:,1),tcav_axis);
    nz(1) = 0;
    nz(end) = 0;
    [ne,eb] = hist(beam(:,2),yag_axis);
    ne(1) = 0;
    ne(end) = 0;
elseif nargin == 6
    [nz,zb] = hist(beam(:,1),Nbin);
    [ne,eb] = hist(beam(:,2),yag_axis);
    ne(1) = 0;
    ne(end) = 0;
else
    [nz,zb] = hist(beam(:,1),Nbin);
    [ne,eb] = hist(beam(:,2),Nbin);
end

ps = hist2(beam(:,1),beam(:,2),zb,eb);

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