function PlotMajorLoops(Mloops,Info)

linecolor = [1 0 0;0 0 1]; % set RGB color to [red;blue]

h = figure;
plot(Mloops.Field{1},Mloops.Kerr_x{1},Mloops.Field{2},Mloops.Kerr_x{2},'Color',linecolor(1,:),'LineWidth',3);
hold all;
plot(Mloops.Field{3},Mloops.Kerr_x{3},Mloops.Field{4},Mloops.Kerr_x{4},'Color',linecolor(2,:),'LineWidth',3);
title('Major Loops comparison');
xlabel('H_{app}');
ylabel('PK');

print(h,'-djpeg','-r300',[Info.DiskID '_MajorLoops']);

%{
delta_mloops.Field{1} = mloops.Field{1};
delta_mloops.Field{2} = mloops.Field{2};
delta_mloops.Kerr_x{1} = mloops.Kerr_x{1}-mloops.Kerr_x{3};
delta_mloops.Kerr_x{2} = mloops.Kerr_x{2}-mloops.Kerr_x{4};

plot(delta_mloops.Field{1},delta_mloops.Kerr_x{1},delta_mloops.Field{2},delta_mloops.Kerr_x{2});
%}

end