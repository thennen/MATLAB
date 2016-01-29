function numout = ndig( numin , ndigits)
% convert positive number to integer string with length ndigits

numin = abs(round(numin));
numout = num2str(numin);

numzeros = ndigits - floor(log10(abs(numin))+1);
numzeros = min(numzeros,ndigits-1); % don't try to add more than ndigits-1 0s

if numzeros < 0;
    error('numin has greater than ndigits')
end

for i=1:numzeros
    numout = ['0' numout];
end

end         