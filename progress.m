
fprintf(1,'Progress:   0/100');
for i=1:100
digits = floor(log10(i)) + 1;
fprintf(1,'\b\b\b\b\b\b\b\b %3d/100',i); pause(.1)
end
fprintf('\n')