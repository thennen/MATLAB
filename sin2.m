function F = sin2(P,xdata)
%function for fitting Asin(2x+B) to some data.  xdata should be given in
%degrees.
A = P(1);
B = 0;
%C = P(3);

F = A*sin(2*(xdata+B)*pi/180);

end

