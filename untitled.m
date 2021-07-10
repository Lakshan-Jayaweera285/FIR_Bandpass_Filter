%------------------------------------------------------
%Defining filter specifications for 140034
%A=0,B=3,C=4

Ap=0.10             %max passband ripple
Aa=53               %max stopband attanuation
omega_p1=700        %Lower passband edge
omega_p2=1000       %upper passband edge
omega_a1=600        %Lower stopband edge
omega_a2=1150       %Upper stopband edge
omega_s=2800        %Sampling frequency

%--------------------------------------------------------
%finding parameters to define window function
%--------------------------------------------------------
%Finding delta and Aa

delta_p=(10.^(0.05*Ap)-1)/(10.^(0.05*Ap)+1)
delta_a=10.^-(0.05*Aa)

delta=min(delta_p,delta_a)

Aa=-20*log10(delta)
%------------------------------------------------------------
%choosing alpha
if Aa<=21
    alpha=0
elseif (Aa<=50)
    alpha=0.5842*((Aa-21).^0.4)+0.07866*(Aa-21)
else
    alpha=0.1102.*(Aa-8.7)
end

%-------------------------------------------------------------
%Choosing parameter D

if Aa<=21
    D=0.9222
else
    D=(Aa-7.95)/14.36
end
%-------------------------------------------------------------
%selecting lowest value for N

Bt=min((omega_p1-omega_a1),(omega_a2-omega_p2))
omega_c1=(omega_p1-Bt/2)
omega_c2=(omega_p2+Bt/2)

N=ceil(omega_s*D/Bt+1)

if rem(N,2)==0                  % if N is an even number,N=N+1 to make it odd
    N=N+1
end

n=-(N-1)/2:(N-1)/2
%--------------------------------------------------------------
%finding KAIZER window function

% window_fn=kaizer_window(N,alpha)

%--------------------------------------------------------------
%defining public variables to check kaizer window
global beta_value
global Io_alpha
global Io_beta

T=2*pi/omega_s
%---------------------------------------------------------------
%defining h(nT)
impulse_res=[]
for n=0:(N-1)/2
    if n==0
        impulse_res(n+1)=2/omega_s*(omega_c2-omega_c1);
    else
        impulse_res(n+1)=1/((n+1)*pi)*(sin(omega_c2*(n+1)*T)-sin(omega_c1*(n+1)*T));
    end
end

impulse_flip=flip(impulse_res)

impulse_response=[impulse_flip(1:length(impulse_flip)-1),impulse_res]

%---------------------------------------------------------------
Window_impulse_fn=impulse_response;
[freq_response,freq]=freqz(Window_impulse_fn);

plot(freq/T,-20*log10(abs(freq_response)));




