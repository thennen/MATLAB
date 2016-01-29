function PlotMinorLoops(MajLoops,MinLoops,Info)

[~,majloopsnum] = size(MajLoops.Kerr_x); 
[~,minloopsnum] = size(MinLoops.Kerr_x); 

% NOT NORMALIZED PLOT
%{
figure;
plot(MajLoops.Field{1},MajLoops.Kerr_x{1},MajLoops.Field{2},MajLoops.Kerr_x{2},'Color','k','LineWidth',2);
hold all;

for i = 1:loopsnum
    plot(MinLoops.Field{i},MinLoops.Kerr_x{i},'LineWidth',2);
    hold all;
end
%}

% NORMALIZED PLOT
% Define and initialize vector to store parameters for normalization
NormFactors = [0 0];

% Using first half major loop, calculates amplitude and offset values needed for
% normalization BY FITTING A HORIZONTAL LINE IN SATURATION REGIME (HIGH FIELD)
%NormFactors = NormalizeRemDecay(MajorLoops);

% Using first half major loop, calculates amplitude and offset values needed for
% normalization BY USING MAX AND MIN PK VALUES
%NormFactors(1) = (max(MajorLoops.Kerr_x{1}) - min(MajorLoops.Kerr_x{2}))/2;     % Amplitude (half peak-to-peak)
%NormFactors(2) = (max(MajorLoops.Kerr_x{1}) + min(MajorLoops.Kerr_x{2}))/2;     % Offset

% Using first half major loop, calculates amplitude and offset values needed for
% normalization BY CONSIDERING VALUES AT H = 0
NormFactors(1) = (MajLoops.Kerr_x{1}(185) - MajLoops.Kerr_x{2}(185))/2;      % Amplitude (half peak-to-peak)
NormFactors(2) = (MajLoops.Kerr_x{1}(185) + MajLoops.Kerr_x{2}(185))/2;      % Offset


for i = 1:majloopsnum
   MajLoops_Norm.Field{i} = MajLoops.Field{i};
   MajLoops_Norm.Kerr_x{i} = (MajLoops.Kerr_x{i} - NormFactors(2)) ./ NormFactors(1); 
end

for i = 1:minloopsnum
   MinLoops_Norm.Field{i} = MinLoops.Field{i};
   MinLoops_Norm.Kerr_x{i} = (MinLoops.Kerr_x{i} - NormFactors(2)) ./ NormFactors(1);
end

h1 = figure;
plot(MajLoops_Norm.Field{1},MajLoops_Norm.Kerr_x{1},MajLoops_Norm.Field{2},MajLoops_Norm.Kerr_x{2},'Color','k','LineWidth',2);
hold all;

for i = 1:minloopsnum
    plot(MinLoops_Norm.Field{i},MinLoops_Norm.Kerr_x{i},'LineWidth',2);
    hold all;
end

plot([0 0],ylim,'k');
hold all;
plot(xlim,[0 0],'k');
title('Minor Loops');
xlabel('H_{app}');
ylabel('PK');
print(h1,'-djpeg','-r300',[Info.DiskID '_MinorLoops']);

end