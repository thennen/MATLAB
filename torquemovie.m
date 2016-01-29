
for k = 1:60
plot(MAngle,16*(sin(MAngle*pi/180).^2) - 600*20/1000*cos(k*360/60*pi/180 - MAngle*pi/180))
hold all
%axis equal
plot(k*360/60,-20:20)
M(k) = getframe;
hold off
end
movie(M,30,5)

