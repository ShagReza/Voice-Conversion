clc; close all; clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total_FFT175=[]
for num_session=5:5
    for num_subject=191:199%wav W564 nadarim
        num_subject
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
            FFT175=[FFT175;(real(frmfft))'];
            m = m + 256 ;  n = n + 1;
            save( ['C:\Users\Admin\BME\project arshad\data\FFT175\FFT_',num2str(num_session),num2str(num_subject)],'FFT175' );
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xx=[];
dd=[];
m=0;
num_session=4;
subfirst=100;
subend=199;
for num_subject=subfirst:subend
    num_subject
    for num_session=4:5
        x=load(['C:\Users\Admin\BME\project arshad\data\FFT175\FFT_',num2str(num_session),num2str(num_subject)]);
        x=x.FFT175;
        xx=[xx;x];
        dd(1+m:length(x)+m,1:(subend-subfirst+1))=0.01;
        dd(1+m:length(x)+m,num_subject-99)=.99; %%cod goyande az 1 be .95
        m=m+length(x);
    end
end
%%%
totallabel=[];
for num_subject=subfirst:subend
    for num_session=4:5
        num_subject
        label=load(['C:\Users\Admin\BME\project arshad\data\LABLES\Zt',num2str(num_session),num2str(num_subject)]);
        label=label.Z;
        totallabel=[totallabel;label];
    end
end
%hazfe sokot:
P=length(dd);
i=1;
while i < P
    if totallabel(i)==40
        %[maxdd,indmax_dd]=max(dd(i,:));
        dd(i,:)=[];
        xx(i,:)=[];
        totallabel(i)=[];
        i=i-1;
        P=P-1;
        P
    end
    i=i+1;
end
XX=xx; DD=dd;
%%%%%%%%%%%%%%%%%%%%

eta=0.5;
max_error=0.01;
alpha=0.7;
p=1;
k=1;
E=0;
S=0;
n_z3=175;
n_z2=175;
n_z=512;
n_y2=128;
n_a=45;
n_b=100;
n_y1=128;
n_x=512;
n_x0=175;
%%%
v0=0.1*rand(n_x0+1,n_x);
dv0=0.001*rand(n_x0+1,n_x);
v1=0.1*rand(n_x+1,n_y1);
dv1=0.001*rand(n_x+1,n_y1);
w1=0.1*rand(n_y1+1,n_a);
dw1=0.001*rand(n_y1+1,n_a);
u1=0.1*rand(n_y1+1,n_b);
du1=0.001*rand(n_y1+1,n_b);
w2=0.1*rand(n_a+1,n_y2);
dw2=0.001*rand(n_a+1,n_y2);
u2=0.1*rand(n_b+1,n_y2);
du2=0.001*rand(n_b+1,n_y2);
v2=0.1*rand(n_y2+1,n_z);
dv2=0.001*rand(n_y2+1,n_z);
v3=0.1*rand(n_z+1,n_z2);
dv3=0.001*rand(n_z+1,n_z2);
v4=0.1*rand(n_z2+1,n_z3);
dv4=0.001*rand(n_z2+1,n_z3);
E=5; k=1; e=[];
%%%
P=length(xx);
eta=(-0.2/500)*k+(0.3+0.2/500);
ESPEACH=[];
ESPEAKER=[];
while E>0.01
    E=0; Espeach=0; Espeaker=0;
    r=randperm(P);
    xx = xx(r,:);
    dd = dd(r,:);
    for p=1:P
        x0=xx(p,:);
        d=dd(p,:);
        x=logsig([x0 1]*v0);
        y1=logsig([x 1]*v1);
        b=logsig([y1 1]*u1);
        a=logsig([y1 1]*w1);
        y2=logsig([a 1]*w2+[b 1]*u2);% d or b!!??
        z=logsig([y2 1]*v2);
        z2=logsig([z 1]*v3);
        z3=[z2 1]*v4;
        %%%
        E=E+sum((x0-z3).^2)+sum((d-b).^2);
        Espeach=sum((x0-z3).^2)+Espeach;
        Espeaker=Espeaker+sum((d-b).^2);
        AA=(+0.5/500)*k+(0.01-0.5/500);
        dz3=(x0-z3)*AA;%*0.01;
        Z2=[z2 1];
        dv4=eta.*Z2'*dz3+alpha.*dv4;
        v4=v4+dv4;
        dZ2=Z2.*(1-Z2).*(dz3*v4');
        dZ2(n_z2+1)=[];
        Z=[z 1];
        dv3=eta.*Z'*dZ2+alpha.*dv3;
        v3=v3+dv3;
        dZ=Z.*(1-Z).*(dZ2*v3');
        dZ(n_z+1)=[];
        Y2=[y2 1];
        dY2=Y2.*(1-Y2).*(dZ*v2');
        dY2(n_y2+1)=[];
        dv2=eta.*Y2'*dZ+alpha.*dv2;
        v2=v2+dv2;
        A=[a 1];
        dA=A.*(1-A).*(dY2*w2');
        dA(n_a+1)=[];
        dw2=eta.*A'*dY2+alpha.*dw2;
        w2=w2+dw2;
        D=[d 1];
        du2=eta.*D'*dY2+alpha.*du2;
        u2=u2+du2;
        db=(d-b).*b.*(1-b);
        Y1=[y1 1];
        dw1=eta.*Y1'*dA+alpha.*dw1;
        w1=w1+dw1;
        du1=eta.*Y1'*db+alpha.*du1;
        u1=u1+du1;
        dY1=Y1.*(1-Y1).*(dA*w1'+db*u1');
        X=[x 1];
        dY1(n_y1+1)=[];
        dv1=eta.*X'*dY1+alpha.*dv1;
        v1=v1+dv1;
        dX=X.*(1-X).*(dY1*v1');
        dX(n_x+1)=[];
        X0=[x0 1];
        dv0=eta.*X0'*dX+alpha.*dv0;
        v0=v0+dv0;
    end

    ESPEACH(k)=sqrt((Espeach)/(length(xx)*18));
    ESPEAKER(k)=sqrt((Espeaker)/(length(xx)*18));
    EE=sqrt(E/(length(xx)*18))
    e(k)=EE;
    k=k+1
    if k<470
        eta=(-0.3/500)*k+(0.3+0.3/500)
    else
        eta=0.01;
    end

    save('C:\Users\Admin\BME\project arshad\data\FFT175\FFt175.mat')

    %hold on;
    %plot(e);
    %drawnow;
end
%%%%%%%%%%%%%%%%


