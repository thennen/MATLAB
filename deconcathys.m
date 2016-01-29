function [Field,Hys] = deconcathys(ConcatField,ConcatHys)
% DECONCATHYS Separates hysteresis loops in a concatenated array of loops


FieldLength = length(ConcatField);
HysLength = length(ConcatHys);

if FieldLength ~= HysLength
    error('Field array has different length than Hysteresis array');
end

j=2;
breakpoint(1) = 1;
for i = 2:FieldLength-1
    % Search for indices where field direction reverses, several options for breakpoint definition
    %if ConcatField(i-1) < ConcatField(i) && ConcatField(i) > ConcatField(i+1);  % increase to decrease
    %if ConcatField(i-1) > ConcatField(i) && ConcatField(i) < ConcatField(i+1);  % decrease to increase
    if (ConcatField(i-1) < ConcatField(i) && ConcatField(i) > ConcatField(i+1)) || (ConcatField(i-1) > ConcatField(i) && ConcatField(i) < ConcatField(i+1));  % any direction change
    %if ConcatField(i-1) < 0 && ConcatField(i) >= 0    % Crosses zero from left to right
    breakpoint(j) = i;
    j = j+1;
    end
end
breakpoint(j) = FieldLength;

for k = 1:j-1
    % Cell Form
    Field{k} = ConcatField(breakpoint(k):breakpoint(k+1));
    Hys{k} = ConcatHys(breakpoint(k):breakpoint(k+1));
    % Matrix Form
    %Field(1:breakpoint(k+1)-breakpoint(k),k) = ConcatField(breakpoint(k):breakpoint(k+1)-1);
    %Hys(1:breakpoint(k+1)-breakpoint(k),k) =  ConcatHys(breakpoint(k):breakpoint(k+1)-1);
    %
end

end

