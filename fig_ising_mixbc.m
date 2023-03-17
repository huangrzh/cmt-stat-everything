function fig = fig_ising_mixbc(Ls, ds, mksize, fname, fontsize, z2)
fig = figure; hold on; box on;
dims0 = 0:6;
degs = [1,1,1,2,2,3,4];
for i = 1:size(ds,1)
    mk = 'ro-';
    if and(mod(i,2)==0,z2==2)
        mk = 'r+-';
    end
    plot(1./log(Ls), ds(i,:), mk, 'MarkerSize',mksize);
end

x1 = min(1./log(Ls));
x2 = max(1./log(Ls));
xl = x1 - (x2-x1)*0.1;
xr = x2 + (x2-x1)*0.03;
for i = 1:length(dims0)
    plot([xl,xr], [dims0(i),dims0(i)], 'k--')
end

xlabel('1/log(\xi)');
fig.OuterPosition(3) = fig.OuterPosition(4);
ax = gca();
ax.YLim(2) = 7.5;
xr2 = x2 + (x2-x1)*0.12;
ax.XLim = [xl,xr2]; 

LA = axis;
text(LA(1)+(LA(2)-LA(1))*0.01,(LA(4)-LA(3))*0.95+LA(3),fname, ...
    'FontSize',fontsize,'FontName','Times New Roman');
ax.FontName = 'Times New Roman';
ax.FontSize = fontsize;


y = cell(7,1);
y{1} = '1/16';
for id = 2:7
    y{id} = [num2str((id-1)*16+1),'/16'];
end
yticks(dims0);
yticklabels(y);

fig_shift_upper(-0.1);
lax = axis;
yyaxis right;
axis(lax);
yticks(dims0+1/16);
y = yticklabels();
x0 = xr;
for i = 1:length(dims0)
    text(x0, dims0(i)+1/16, ['x',num2str(z2*degs(i))], 'Color', 'r',...
        'FontSize', fontsize-4);
    y{i} = '[\sigma]';
end
yticklabels(y);
ax.YColor = 'r';
ax = ax.YAxis(2);
ax.TickDirection = 'none';
end

function fig_shift_upper(x)
lax = axis;
lax(3) = x;
axis(lax);
end