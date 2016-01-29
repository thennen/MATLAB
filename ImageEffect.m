function out = ImageEffect( in , range)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Field = in.Field;
X = in.X;
out = in;
fit = polyfit(Field(Field >= range(1) & Field <= range(2)),X(Field >= range(1) & Field <= range(2)),1);
out.Fitline = fit(2) + fit(1) .* Field;
out.Correction = X ./ out.Fitline;


end