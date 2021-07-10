omega = -pi:0.02:pi;
cut_off = pi/2;
syms n
f(n) = (sin(n*cut_off)/(n*pi))*exp(-1i*omega*n);
H_omega = 0;
for i=-100:100
    if i~=0
        H_omega = H_omega + f(i);
    else
        H_omega = H_omega + (cut_off/pi);
    end
end
p = plot(omega,H_omega);
title('Frequency Response for N=100')
xlabel('Omega')
ylabel('H(w)')
grid on
set(p,'Color','red','LineWidth',2)