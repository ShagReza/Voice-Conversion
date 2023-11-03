%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% project arshad %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% 5 layer NN %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc,close all,clear all
% ------loading data for training----------------------
xx=[]; m=0;
num_session=4;
for num_subject=107:10
    x=load(['C:\Users\Admin\BME\project arshad\data\LHCB_ normalized\LHCB_',num2str(num_session),num2str(num_subject)]);
    x=x.CBw;
    xx=[xx;x];
    m=m+length(x);
end
xx(:,19:end)=[];

xx=[]; dd=[];
xx=load(['C:\Users\Admin\BME\project arshad\data\LHCB_ normalized\LHCB_',num2str(num_session),num2str(num_subject)]);
xx=xx.CBw; xx(:,19:175)=[];

%%%
input=xx;
size_select = size(input);
P=size_select(1);
%%%
e=[];
eta=.5;
MaxE=.01;
alpha=.7;
M=18;   %first layer
N=128;    %second layer
K=1;    %middle layer
%%%
Wa=0.01*rand(M+1,N)-0.005;
Va=0.01*rand(N+1,K)-0.005;

Wb=0.01*rand(K+1,N)-0.005;
Vb=0.01*rand(N+1,M)-0.005;

DeltaWa=zeros(M+1,N);
DeltaVa=zeros(N+1,K);

DeltaWb=zeros(K+1,N);
DeltaVb=zeros(N+1,M);

clear CB r size_select
p=1;k=1;E=0;counter=1;F=1;q=1;
tic
%%%%
eta=(-0.4/4999)*k+(0.5+0.4/4999);
while (F) && (eta>0)
    k
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
    eta=(-0.4/4999)*k+(0.5+0.4/4999);
end


% % ------loading data for test----------------------
xx=[]; m=0;
num_session=4;
num_subject=107;
x=load(['C:\Users\Admin\BME\project arshad\data\LHCB_ normalized\LHCB_',num2str(num_session),num2str(num_subject)]);
xx=x.CBw;
xx(:,19:end)=[];
%%%
input=xx;
size_select = size(input);
P=size_select(1);
p=1;
XB=[]; Z=[];
while p<=P
    xa=input(p,:);
    xa=[xa  1];
    ya=logsig(xa*Wa);
    YA(p,:)=ya;
    ya = [ya , +1];
    z=logsig(ya*Va);
    Z(p,:)=z;
    z=[z,1];
    yb=logsig(z*Wb);
    yb = [yb , +1];
    xb=yb*Vb;
    XB(p,:)=xb;
    p=p+1;
end
%%%%%%%%%% jobrane miyangingiri:
% mean_CB=load('Z:\project arshad\programs\p5_mean & var extraction\mean_CB.mat');
% mean_CB=mean_CB.mean_CB;
% var_CB=load('Z:\project arshad\programs\p5_mean & var extraction\var_CB.mat');
% var_CB=var_CB.var_CB;
% 
% mean_CB1=mean_CB(num_subject,:);
% var_CB1=var_CB(num_subject,:);
% 
% XB(:,1:18) = XB(:,1:18) .* repmat(var_CB1, size(XB , 1) ,1) ;
% XB(:,1:18) = XB(:,1:18) + repmat(mean_CB1 , size(XB , 1),1) ;
%%%%%%%%%% bazsazi soot az khoroji shabake %%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
CBFBS2 = CBFBS;
CBFBS2(1:175,19:175)=rand(175,157);
%%%%%%%%%%%%%
CBw=XB;
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









plot(X100(:,3),X100(:,4),'.b'),hold on,plot(X101(:,3),X101(:,4),'.r'),






