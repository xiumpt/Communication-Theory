% =========== comm HW4 ===========
% Q:
% bits: 11001000 (b0b1...b7) using OFDM(fc=20Hz) to transmit
% b0b1 with QPSK modulation (fc=18Hz)
% b2b3 with QPSK modulation (fc=19Hz)
% b4b5 with QPSK modulation (fc=21Hz)
% b6b7 with QPSK modulation (fc=22Hz)
% ================================

%%	(a) plot OFDM signal s(t), 0<=t<=1sec
close all;clear all;
% signal: 11(+1,+1)  00(-1,-1)  10(+1,-1)  00(-1,-1)  


fs=256; % sample frequencys
t=0:1/fs:1; % sample time
fc_s1=18; fc_s2=19; fc_s3=21; fc_s4=22;
A1=1; B1=1;  A2=-1; B2=-1;    A3=1; B3=-1;     A4=-1; B4=-1;

s1=A1*cos(2*pi*fc_s1*t)-B1*sin(2*pi*fc_s1*t);	
s2=A2*cos(2*pi*fc_s2*t)-B2*sin(2*pi*fc_s2*t);	
s3=A3*cos(2*pi*fc_s3*t)-B3*sin(2*pi*fc_s3*t);	
s4=A4*cos(2*pi*fc_s4*t)-B4*sin(2*pi*fc_s4*t);	
 s = s1+s2+s3+s4;
 %{ 
subplot(5,1,1);plot(t,s1);
    title('fc=18Hz');
subplot(5,1,2);plot(t,s2);
    title('fc=19Hz');
subplot(5,1,3);plot(t,s3);
    title('fc=21Hz');
subplot(5,1,4);plot(t,s4);
    title('fc=22Hz');

subplot(5,1,5);plot(t,s);
        title('OFDM signal s(t)');
%}        
figure; plot(t,s);
    title('OFDM signal s(t)');
    xlabel('time'); ylabel('amplitude');
    axis([-inf,inf,-inf,inf]);

%%  (b) plot OFDM signal spectrum S(f)
%{
N=length(t);    %smaple length
df=fs/(N-1);    % frequency resolution
ofdm_f=(0:N-1)*df;  % frequency of every point
Y=fft(s(1:N))/N*2;  % real amplitude

subplot(2,1,1); plot(ofdm_f(1:N/2),abs(Y(1:N/2)));  % single side band
    title('OFDM single side band frequency spectrum S(f)');
    xlabel('frequency'); ylabel('amplitude');
    axis([-inf,inf,-inf,inf]);

subplot(2,1,2); plot(ofdm_f,abs(Y));    % double side band
    title('OFDM double side band frequency spectrum S(f)');
    xlabel('frequency'); ylabel('amplitude');
    axis([-inf,inf,-inf,inf]);
%}
%%  (c) plot OFDM complex envelope spectrum S(f)
%{
[up,down] = envelope(s);
figure;
plot(t,s);hold on;
plot(t,up,t,down,'linewidth',1.5); 
    title('OFDM complex envelope');
    xlabel('time'); ylabel('amplitude');
    
up_fft = fft(up(1:N))/N*2;
up_fftshift = fftshift(up_fft);

figure; plot(ofdm_f,abs(up_fftshift));
    title('OFDM complex envelope spectrum S(f)');
    xlabel('frequency'); ylabel('amplitude');
    axis([-inf,inf,-inf,inf]);
%}
%%  (d) plot the in-phase components s_I(t), 0<=t<=1sec
    %   (e) plot the quadrature component s_Q(t), 0<=t<=1sec
    %   (h) describe how to build the OFDM signal s(t) by  s_I(t) and  s_Q(t),
    %          and plot to compare the signal is same as s(t).
fc_ofdm=20;
figure; 

%   in-phase components 
s_i=2*sin(4*pi*t);
subplot(4,1,1);plot(t,s_i);title('in-phase components s_I(t)');
    xlabel('time'); ylabel('amplitude');
    axis([-inf,inf,-inf,inf]);
    
%   quadrature component
s_q=-2*sin(4*pi*t)+2*sin(2*pi*t)-2*cos(2*pi*t);
subplot(4,1,2);plot(t,s_q);title('quadrature component s_Q(t)');
    xlabel('time'); ylabel('amplitude');
    axis([-inf,inf,-inf,inf]);
    
%   s(t) = s_i*cos + s_q*sin
sx=s_i.*cos(2*pi*fc_ofdm*t)-s_q.*sin(2*pi*fc_ofdm*t);
subplot(4,1,3);plot(t,sx);title('s_I(t)*cos(2\pif_ct) - s_Q(t)*sin(2\pif_ct)');
    xlabel('time'); ylabel('amplitude');
    axis([-inf,inf,-inf,inf]);
    
%   original signal
subplot(4,1,4);plot(t,s);title('OFDM signal s(t)');
    xlabel('time'); ylabel('amplitude');
    axis([-inf,inf,-inf,inf]);
    
%%  (f) plot envelope sqrt(s_I(t)^2+s_Q(t)^2) and OFDM signal s(t) on the same figure to compare
envelope = sqrt(s_i.^2+s_q.^2);

figure;
plot(t,envelope,'linewidth',1.5);hold on;
plot(t,-envelope,'linewidth',1.5);hold on;
plot(t,s);
    text(0.45, 4, '\bf\leftarrow upper envelope = \surd(s_I(t)^2+s_Q(t)^2)');
    text(0.03, -4, '\bf lower envelope\rightarrow');
    text(0.6,0,'\bf\leftarrow OFDM signal');
    title('OFDM envelope & OFDM signal');
    xlabel('time'); ylabel('amplitude');
    axis([-inf,inf,-inf,inf]);
    
%%  (g) describe how to get the sample of s_I(t) and s_Q(t) by 8IFFT,
%            and plot the corresponding envelope and its sample on the same figure. 
t_8ifft=0:1/8:1;

%   in-phase components 
s_i2=2*sin(4*pi*t_8ifft);
subplot(2,1,1);stem(t_8ifft,s_i2);
    title('s_I(t) sample by 8IFFT');
    xlabel('time'); ylabel('amplitude');hold on;
    axis([-inf,inf,-inf,inf]);
    
%   quadrature component
s_q2=-2*sin(4*pi*t_8ifft)+2*sin(2*pi*t_8ifft)-2*cos(2*pi*t_8ifft);
subplot(2,1,2);stem(t_8ifft,s_q2);
    title('s_Q(t) sample by 8IFFT');
    xlabel('time'); ylabel('amplitude');
    axis([-inf,inf,-inf,inf]);
    
%   s(t) = s_i*cos + s_q*sin
sx2=s_i2.*cos(2*pi*fc_ofdm*t_8ifft)-s_q2.*sin(2*pi*fc_ofdm*t_8ifft);
figure; stem(t_8ifft,sx2);    
    hold on; plot(t,sx);
    title('envelope sample by 8IFFT & envelope signal');
    xlabel('time'); ylabel('amplitude');
    axis([-inf,inf,-inf,inf]);    
    
