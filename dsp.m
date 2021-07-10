clear all;
close all;
w = -2.5*pi:0.01:2.5*pi;
y = exp(2i*w);
temp = sin((2.5)*w)./sin((0.5)*w);
x = y.*temp;
x1 = abs(x);
plot(w,x1);grid;
title('dgss');


