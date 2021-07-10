%%
close all;
clear all;
clc

%%this is the kaiser window method for band pass filter
% my index number is 180288L

% getting A=2, B=8 and C=8 
A=2;
B=8;
C=8;



%Basic parameters for the Band Pass filter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_p = 0.03+0.01*A;  %% (dB)max passband ripple
A_a = 45+B;          %%(dB) min stopband attenuation
omega_p1 = C*100+300;     %% rad/s  , lower passband edge
omega_p2 = C*100+700;     %%rad/s  , upper passband edge
omega_a1 = C*100+150;     %%rad/s , lower stopband edge
omega_a2 = C*100+800;     %%rad/s , upper stopband edge 
omega_s = 2*(C*100+1200);  %%rad/s  , sampling frequency
fs = omega_s/(2*pi);      %(Hz)  sampling freq
Len = 275;% no of samples for plot input signals(in the assingment noted to get between 200 to 300 value for this)
T = 1/fs;    %%period of sampling


%%%%%%%%%%%%%%%%%%%%%%%%%%Kaiser Window Generation%%%%%%%%%%%%%%%%%%%%%%%
%%
% get critical transition width
Bt = min((omega_p1-omega_a1),(omega_a2-omega_p2));

omega_c1=omega_p1-Bt/2;   %rad/s ,this is the ideal freq of lower cutoff 
omega_c2= omega_p2+Bt/2;  %rad/s ,this is the ideal freq of upper cutoff 

%actual passband ripple
delta_P=(10^(0.05*A_p)-1)/(10^(0.05*A_p)+1);
delta_A=10^(-0.05*A_a);

delta=min(delta_P,delta_A); %min value we get to delta value in next steps

%actual attenuation using delta
Aa = - 20*log10(delta); %(dB)

%get parameter alpha value using Aa 
%these are the defined values so it has three main boundaries
if Aa<=21
    alpha=0;
elseif Aa<=50
    alpha=0.5842*(Aa-21)^0.4+0.07886*(Aa-21);
else
    alpha=0.1102*(Aa-8.7);
end

%get parameter D using Aa value
if Aa<=21
    D=0.9222;
else
    D=(Aa-7.95)/14.36;
end

%select N- this is the temperary value , it makes N value
temp_var=(omega_s*D/Bt)+1;
if (temp_var-fix(temp_var)) ==0
    if mod(temp_var,2)==1
        N=temp_var;
    else
        N=temp_var+1;
    end
else
    if mod(fix(temp_var),2)==0
        N=fix(temp_var)+1;
    else
        N=fix(temp_var)+2;
    end
end

%get M
M=(N-1)/2;     %this M is value for one side number of discrete signals

%set non zero domain for window 
n=linspace(-M,M,2*M+1);

%get the beta values
beta=alpha*sqrt(1-(n/M).^2);

numerator=getBessel(beta);  % return bessel function value using the beeta
denominator=getBessel(alpha); % return bessel function value using alpha

%kiaser window generation
w_n=numerator/denominator;  % equation of the kaiser window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%generate ideal bandpass filter
h_n=sin(omega_c2*n/fs)./(pi*n)-sin(omega_c1*n/fs)./(pi*n);
h_n(M+1)=(omega_c2-omega_c1)*2/omega_s;

%FIR BP Filter
hw_n=h_n.*w_n;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%calculations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
[h,f]=freqz(hw_n);%frequency response of filter in normalized freq range

f=f*(omega_s/(2*pi));%freq domain correction
H_db=20*log10(abs(h));%get the magnitude in db scale

%causal impulse response of the filter
figure;
stem(n(M+1:2*M+1),hw_n(M+1:2*M+1));
xlabel('n (Samples)');ylabel('Amplitude');title('impulse response of the filter')

minPassRip=-A_p/2*ones(1,length(f));
maxStopRip=-A_a*ones(1,length(f));
maxPassRip=A_p/2*ones(1,length(f));

%frequency response of the filter
figure;
plot(f,H_db);
hold on;
plot(f,minPassRip);
hold on;
plot(f,maxPassRip);
hold on;
plot(f,maxStopRip);
xlabel('Frequency (rad/s)');ylabel('Magnitude (dB)');title('bandpass FIR filter (Kaiser Windowed)frequency response')




figure;
plot(f,abs(h));
xlabel('Frequency (rad/s)');ylabel('Magnitude');title('BandPass filter');
%%%%%%%%%%%%%%%%%%%%%%validating results%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
n=(0:Len-1);%integer domain
f=(-Len/2:Len/2-1)*(omega_s/Len);%frequency domain

%input signal
x_n=sin(n*(omega_c1+omega_c2)/(2*fs))+sin(n*(0+omega_a1)/(2*fs))+sin(n*(omega_a2+omega_s/2)/(2*fs));
t_n=sin(n*(omega_c1+omega_c2)/(2*fs));%ideal output

%%%%%%%%%%%%%%%%%%%%%%%frequency domain calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%generate discrete fourier transform of the signals
X_n=fftshift(fft(x_n,Len));     % DFT to x_n using fast fourier transform
H_n=fftshift(fft(hw_n,Len));    % DFT to h_n using fast fourier transform
T_n=fftshift(fft(t_n,Len));     % DFT to t_n using fast fourier transform

%convolute time domain signals
Y_n=(X_n.*H_n);   % in frequency domain multiplication

%plot the graphs in frequency domain 
figure;
subplot(3,1,1);plot(f,abs(X_n/Len));
xlabel('Frequency in (rad/s)');ylabel('Magnitude');title('Input signal');
subplot(3,1,2);plot(f,abs(H_n));
xlabel('Frequency in (rad/s)');ylabel('Magnitude');title('Band pass filter');
subplot(3,1,3);plot(f,abs(Y_n/Len));
xlabel('Frequency in(rad/s)');ylabel('Magnitude');title('Output signal');

figure;
subplot(2,1,1);plot(f,abs(Y_n/Len));
xlabel('Frequency (rad/s)');ylabel('Magnitude');title('Output of the FIR BPF(Frequency Domain)');
subplot(2,1,2);plot(f,abs(T_n/Len));
xlabel('Frequency (rad/s)');ylabel('Magnitude');title('Output of the Ideal BPF(Frequency Domain)');




%%%%%%%%%%%%%%%%%%%%%%%%%generate time domain results and plots%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%arrange the filter
hw_n1=zeros(1,Len);
hw_n1(floor(Len/2)-M:floor(Len/2)+M)=hw_n;

%filter the signal
y_n=conv(x_n,hw_n,'full');

%plot the graphs in time domain
figure;
subplot(3,1,1);plot(n,x_n);
ylim([min(x_n) max(x_n)]);xlim([0 Len]);xlabel('Samples (n)');ylabel('Amplitude');title('Input signal');
subplot(3,1,2);plot(n,hw_n1);
ylim([min(hw_n1) max(hw_n1)]);xlim([0 Len]);xlabel('Samples (n)');ylabel('Amplitude');title('FIR BPF(time domain) ');
subplot(3,1,3);plot(n,y_n(M:Len+M-1));
ylim([-1 1]);xlim([0 Len]);xlabel('Samples (n)');ylabel('Amplitude');title('FIR BPF Output');

figure;
subplot(2,1,1);plot(n,y_n(M:Len+M-1)); %BPF
xlabel('Samples (n)');ylabel('Amplitude');title('FIR BPF Output');ylim([-1 1]);xlim([0 Len]);
subplot(2,1,2);
plot(n,t_n);
ylim([-1 1]);xlim([0 Len]);xlabel('Samples (n)');ylabel('Amplitude');title('ideal BPF Output');

