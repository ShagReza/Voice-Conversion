%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% wave to FFT(175)_without normalization &  fs=22050 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;  clear all; close all
total_FFT175=[]
for num_session=4:5
    for num_subject=1:1 %wav W564 nadarim
        filename = ['C:\Users\Admin\BME\project arshad\data\WAV\W',num2str(num_session),num2str(num_subject),'.wav'];
        [s1,fs] = wavread(filename);
        %%%
        %-feature extraction--------------------------------------------------
        N = 512 ; %1024		            % Frame length
        stp = 256;	%512	            % Frame rate
        L = 18 ;		            % Number of LHCB parameters
        Q = 18 ;		            % Number of Filters
        fc = 7500 ;	                % Final frequency is 7500 Hz
        fs = 22050 ;  %44100          	% Sampling frequency
        kc = round( fc * N / fs );	% Equivallent of 7500 Hz in DFT domain
        f = fs / N * (0:kc);		% Frequencies samples in hertz
        win =[0.54+0.46*cos(2*pi/511*(-256:255))]';	   % Temporal window
        fbark = 6 * log( f / 600 + sqrt( ( f / 600 ).^2 + 1 ) ) ;
        CBFBS1 = zeros(Q,length(f));

        for zk = 1 : Q
            bk = ( ( fbark >= zk - 1 ) & ( fbark < zk +1 ) ) ;
            CBFBS1(zk , :) = bk .* ( 0.5 + 0.5 * cos( pi * ( fbark - zk ) ) ) ;
        end

        CBFBS1 = CBFBS1'; CBFBS2 = CBFBS1 .^ 2 ;
        CBFBS2(1:175,19:175)=rand(175,157); CBFBS=CBFBS2;
        CBw = []; CB2 = [];
        m = 1; n=1; x=1;
        cbfinal1=[]; cbfinal2=[];
        FFT175=[];
        while m < length(s1) - N + 1,
            frm = s1 ( m : m + N - 1 );
            frm = frm - mean(frm);
            frmwin = frm .* win;
            frmfft = fft ( frmwin , N ) ;
            frmfft ( kc+2 : end ) = [];
            frmfftsqr =(( real ( frmfft .* conj ( frmfft ) ))*4)';
            %energy of filter
            %cb = frmfftsqr * CBFBS ;
            %cbfinal1=[cbfinal1 cb];
            %cbcepst = log ( 0.01 + cb ) ;
            %cbfinal2=[cbfinal2 cbcepst];
            %after log
            %CBw = [CBw ; cbcepst] ;
            FFT175=[FFT175;(real(frmfft))'];
            m = m + 256 ;  n = n + 1;
        end
        total_FFT175=[total_FFT175;FFT175];
        %save( ['C:\Users\Admin\BME\project arshad\data\LHCB_without normalization\LHCB_',num2str(num_session),num2str(num_subject)],'CBw' );
        clear CB;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%
input=total_FFT175;
size_input = size(input);
P=size_input(1);
%%%
e=[];
eta=.5;
MaxE=.01;
alpha=.7;
M=175;   %first layer
N=256;    %second layer
K=50;    %middle layer
%%%
Wa=0.01*rand(M+1,N)-0.005;
Va=0.01*rand(N+1,K)-0.005;

Wb=0.01*rand(K+1,N)-0.005;
Vb=0.01*rand(N+1,M)-0.005;

DeltaWa=zeros(M+1,N);
DeltaVa=zeros(N+1,K);

DeltaWb=zeros(K+1,N);
DeltaVb=zeros(N+1,M);

clear CB r size_input
p=1;k=1;E=0;counter=1;F=1;q=1;
tic
%%%%
eta=(-0.4/3000)*k+(0.5+0.4/3000);
while (F)
    r=randperm(P);
    input = input(r,:);
    d = input;   %desired output
    while p<=P
        xa=input(p,:);
        xa=[xa  1];
        ya=logsig(xa*Wa);
        ya = [ya , +1];
        z=logsig(ya*Va);
        z = [z , +1];
        yb=logsig(z*Wb);
        yb = [yb , +1];
        xb=(yb*Vb);

        E=E+sum((d(p,:)-xb).^2);

        dxb=(d(p,:)-xb)*.001;
        dyb=(dxb*Vb') .* yb .* (1-yb);
        dyb(end) = [];
        dz=(dyb*Wb') .* z .* (1-z);
        dz(end) = [];
        dya=(dz*Va') .* ya .* (1-ya);
        dya(end) = [];

        DeltaVb=eta.* yb' *dxb+alpha.*DeltaVb;
        Vb=Vb+DeltaVb;
        DeltaWb=eta.* z' *dyb+alpha.*DeltaWb;
        Wb=Wb+DeltaWb;
        DeltaVa=eta.* ya' *dz+alpha.*DeltaVa;
        Va=Va+DeltaVa;
        DeltaWa=eta.* xa' *dya+alpha.*DeltaWa;
        Wa=Wa+DeltaWa;
        p=p+1;

    end
    e(k)=sqrt(E/(M*P));
    k
    H=e(k)
%     hold on;
%     plot(e);
%     drawnow;
    if H <=MaxE
        F=0;
    else
        E=0;
        k=k+1;
        p=1;
    end
    eta=(-0.4/3000)*k+(0.5+0.4/3000)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=1;
P=length(FFT175);
 while p<=P
        xa=FFT175(p,:);
        xa=[xa  1];
        ya=logsig(xa*Wa);
        ya = [ya , +1];
        z=logsig(ya*Va);
        z = [z , +1];
        yb=logsig(z*Wb);
        yb = [yb , +1];
        xb=(yb*Vb);
        X(p,:)=xb;
        p=p+1;
 end
 %%
num_session=5;
num_subject=1;
%pitch:
filename = ['C:\Users\Admin\BME\project arshad\data\WAV\W',num2str(num_session),num2str(num_subject),'.wav'];
[s1,fs] = wavread(filename);
m = 1; PP=0; PITCH=[]; 
while m < length(s1) - 512 + 1,
    frm = s1( m : m + 512 - 1 );
    a = real ( ifft ( abs ( fft ( frm , 512 ) ) ) );
    [cc,bb] = max ( a ( 50 : 256 ) ) ; %(100:512)
    pitch = bb + 49 ; %99
    m = m + 256; %m=m+512
    PITCH=[PITCH pitch];
end
mm=0;pitch_new=[];jj=0;
for i=1:length(PITCH)
    while mm<256 %512
        mm=mm+PITCH(i);
        if mm<256 %512
            jj=jj+1;
            pitch_new(jj)=PITCH(i);
        else
            pitch_new(jj+1)=256; %512
            pitch_new(jj+2)=mm-256; %512
        end
    end
    mm=mm-256;%512
    jj=jj+2;
end
Pitch=pitch_new;
%%%%%%%%%
fs = 22050 ; fc = 7500 ;	
kc = round( fc * 512 / fs );	
f = fs / 512 * (0:kc);
fbark = 6 * log( f / 600 + sqrt( ( f / 600 ).^2 + 1 ) ) ;
CBFBS1 = zeros(18,length(f));
for zk = 1 : 18
    bk = ( ( fbark >= zk - 1 ) & ( fbark < zk +1 ) ) ;
    CBFBS1(zk , :) = bk .* ( 0.5 + 0.5 * cos( pi * ( fbark - zk ) ) ) ;
end
CBFBS1 = CBFBS1'; CBFBS = CBFBS1 .^ 2 ;
CBFBS2 = CBFBS1;
CBFBS2(1:175,19:175)=rand(175,157);
%%%%%%%%%
frmfftsqr =(( real ( X' .* conj ( X' ) ))*4)';
cb = frmfftsqr * CBFBS ;
cbfinal1=[cbfinal1 cb];
cbcepst = log ( 0.01 + cb ) ;
cbfinal2=[cbfinal2 cbcepst];
CBw = [CBw ; cbcepst] ;
%%%%%%%%%%%%%%%%%%%%%%%%%%
CBw=X;
m1=1; N2=175; m2=1;
newfrms = exp ( CBw ) - 0.01 ;
z=size(newfrms,1)*512-(size(newfrms,1)-1)*256;
smp=zeros(1,z);
PP=0;  q=0;
%kk = newfrms( m2,:) / CBFBS2; %!!!
%kk = newfrms( m2,:) / CBFBS; %new_khodam
kk = newfrms( m2,:)*pinv(CBFBS);  %/(CBFBS);
kk ( 512 : -1 : 339 ) = kk ( 2 : 175 );  % mirrorong
kk ( 176 : 340 ) = 0;
dd =  real ( ifft ( sqrt(kk) , 512 ) ) ;
d = dd .* hamming(512)' ;  %for deletting the effect of needle picks at first and end of frames
d2 = [d(257:512),d(1:256)];
smp ( 256 * q + 1 : 256 * q + 512 ) = smp ( 256 * q + 1 : 256 * q + 512 ) + d2 ;

while ( Pitch ( m1 ) ~= 256 ) & ( m2 == 1 )
    m1=m1+1;
end

m2 = m2 + 1; m1 = m1 + 1;

while ( m2 <= size(newfrms,1) ) & ( m1 < size ( Pitch , 2 ) )
    %kk = newfrms( m2,:) / CBFBS2; %!!!
    %kk = newfrms( m2,:) / CBFBS; %new_khodam
    kk = newfrms( m2,:)*pinv(CBFBS);  %/(CBFBS);
    kk ( 512 : -1 : 339 ) = kk ( 2 : 175 );  % mirrorong
    kk ( 176 : 340 ) = 0;
    dd =  real ( ifft ( sqrt(kk) , 512 ) ) ;
    d = dd .* hamming(512)' ;  %for deletting the effect of needle picks at first and end of frames
    d2 = [d(257:512),d(1:256)];
    PP = Pitch ( m1 );
    while ( Pitch ( m1 ) ~= 256 ) & ( m1 < size ( Pitch , 2 ) )
        smp ( 256 * q + PP : 256 * q + 512 + PP  -1) = smp ( 256 * q + PP : 256 * q + 512 + PP  -1) + d2 ;
        m1 = m1 + 1 ;
        PP = PP + Pitch ( m1 ) -1 ;
    end
    q = q + 1;
    m2 = m2 + 1;
    m1 = m1 + 1;
    newfrms11=newfrms(:,1:18).^2;
    newfrms1=sum( newfrms11);
    d211=d2.^2;
    d21=sum(d211);
end


input.layer=18
1st.layer=128
2nd.layer=(32-5)
3th.layer=128
4th.layer=18
output.layer=18


