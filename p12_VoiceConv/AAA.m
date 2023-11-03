%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% In the name of God %%%%%%%
%%%%%%% shaghayegh reza %%%%%%%%%
%%%%%%%% project arshad %%%%%%%%%
%%%%% 18 khordad 87 %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear all;
%dadegan normalize
%code goyande: .95 & .05
%hazfe sokot:label 40
%shabake 8 laye
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xx=[];
dd=[];
m=0;
num_session=4;
subfirst=100;
subend=199;
for num_subject=subfirst:subend
    for num_session=4:5
        x=load(['C:\Users\Admin\BME\project arshad\data\LHCB_ normalized\LHCB_',num2str(num_session),num2str(num_subject)]);
        x=x.CBw;
        xx=[xx;x];
        dd(1+m:length(x)+m,1:(subend-subfirst+1))=0.01;
        dd(1+m:length(x)+m,num_subject-99)=.99; %%cod goyande az 1 be .95
        m=m+length(x);
    end
end
xx(:,19:end)=[];
%%%
totallabel=[];
for num_subject=subfirst:subend
    for num_session=4:5
        label=load(['C:\Users\Admin\BME\project arshad\data\LABLES\Zt',num2str(num_session),num2str(num_subject)]);
        label=label.Z;
        totallabel=[totallabel;label];
    end
end
%just vowel AA
P=length(dd); j=0; DD=[]; XX=[];
for i=1:P
    if totallabel(i)==1
        j=j+1;
        DD(j,:)=dd(i,:);
        XX(j,:)=xx(i,:);
    end
end
xx=XX; dd=DD;
%%%%%%%%%%%%%%%%%%%%
eta=0.5;
max_error=0.01;
alpha=0.7;
p=1;
k=1;
E=0;
S=0;
n_z3=18;
n_z2=18;
n_z=200;
n_y2=120;
n_a=45;
n_b=100;
n_y1=120;
n_x=200;
n_x0=18;
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
        y2=logsig([a 1]*w2+[d 1]*u2);% d or b!!??
        z=logsig([y2 1]*v2);
        z2=logsig([z 1]*v3);
        z3=[z2 1]*v4;
        %%%
        E=E+sum((x0-z3).^2)+sum((d-b).^2);
        Espeach=sum((x0-z3).^2)+Espeach;
        Espeaker=Espeaker+sum((d-b).^2);
        %AA=(+0.5/500)*k+(0.01-0.5/500);
        AA=0.01;
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
        %db=(d-b).*b.*(1-b);
        db=(d-b)*0.01;
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

    save('C:\Users\Admin\BME\project arshad\programs\p12_VoiceConv\AA.mat')

    %hold on;
    %plot(e);
    %drawnow;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_subject=122;
Esubject=[];
%for num_subject=100:199
xxt=[];
ddt=[];
Esub=0;
m=0;
num_session=4;
xxt=load(['C:\Users\Admin\BME\project arshad\data\LHCB_ normalized\LHCB_',num2str(num_session),num2str(num_subject)]);
xxt=xxt.CBw;
ddt(1+m:length(xxt)+m,1:(subend-subfirst+1))=0.01;
ddt(1+m:length(xxt)+m,num_subject-99)=.99; %%cod goyande az 1 be .95
xxt(:,19:end)=[];
%%%%%%%%%% bazsazi soot az khoroji shabake %%%%%%%%%%
B=[]; A=[];
P=length(xxt);
for p=1:P
    x0=xxt(p,:);
    d=ddt(p,:);
    x=logsig([x0 1]*v0);
    y1=logsig([x 1]*v1);
    b=logsig([y1 1]*u1);
    B(p,:)=b;
    a=logsig([y1 1]*w1);
    y2=logsig([a 1]*w2+[b 1]*u2);
    z=logsig([y2 1]*v2);
    z2=logsig([z 1]*v3);
    z3=[z2 1]*v4;
    Z3(p,:)=z3;
    %Esub=Esub+sum((z3-x0).^2)/(P*18);
    %Esubject(num_subject)=Esub; %mohasebeye khata baraye har goyande be tore mojaza
end
%end
CBw(:,1:18)=Z3(:,1:18);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_CB=load('C:\Users\Admin\BME\project arshad\programs\p5_mean & var extraction\mean_CB.mat');
mean_CB=mean_CB.mean_CB;
var_CB=load('C:\Users\Admin\BME\project arshad\programs\p5_mean & var extraction\var_CB.mat');
var_CB=var_CB.var_CB;

mean_CB1=mean_CB(num_subject,:);
var_CB1=var_CB(num_subject,:);

CBw(:,1:18) = CBw(:,1:18) .* repmat(var_CB1, size(CBw , 1) ,1) ;
CBw(:,1:18) = CBw(:,1:18) + repmat(mean_CB1 , size(CBw , 1),1) ;
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%be dast avardane mizane sehat tashkhise goyande baraye avahaye mokhtalef:
totallabel=[];
for num_subject=100:104
    label=load(['C:\Users\Admin\BME\project arshad\data\LABLES\Zt',num2str(num_session),num2str(num_subject)]);
    label=label.Z;
    totallabel=[totallabel;label];
end

accuracy(1,43)=0;
num(1,43)=0;
for i=1:P
    [maxB,indmax_B]=max(B(i,:)) %indmax:indice of max B
    [maxdd,indmax_dd]=max(dd(i,:));
    LABEL=totallabel(i);
    num(LABEL)=num(LABEL)+1;
    if indmax_B==indmax_dd
        accuracy(LABEL)=accuracy(LABEL)+1;
    end
end
percent_accuracy=[];
for i=1:43
    percent_accuracy(i)=accuracy(i)/num(i)*100
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





input.layer=18
1st.layer=100
2nd.layer=60
3th.layer=(32-5)
4th.layer=60
5th.layer=100
6th.layer=60
output.layer=18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Espeaker(1:45)=0; Ephonem(1:45)=0; num(1:45)=1;  subiden(1:45)=0; speakeriden(1:199)=0;
for num_subject=100:199
    L=0;
    for num_session=4:5
        xxt=[];
        ddt=[];
        Esub=0;
        m=0;
        num_session=4;
        xxt=load(['C:\Users\Admin\BME\project arshad\data\LHCB_ normalized\LHCB_',num2str(num_session),num2str(num_subject)]);
        xxt=xxt.CBw;
        ddt(1+m:length(xxt)+m,1:(subend-subfirst+1))=0.01;
        ddt(1+m:length(xxt)+m,num_subject-99)=.99; %%cod goyande az 1 be .95
        xxt(:,19:end)=[];
        label=load(['C:\Users\Admin\BME\project arshad\data\LABLES\Zt',num2str(num_session),num2str(num_subject)]);
        label=label.Z;
        B=[]; A=[];
        P=length(xxt);
        L=L+P;
        for p=1:P
            x0=xxt(p,:);
            d=ddt(p,:);
            x=logsig([x0 1]*v0);
            y1=logsig([x 1]*v1);
            b=logsig([y1 1]*u1);
            B(p,:)=b;
            [a1,b1]=max(b);
            if b1==(num_subject-99)
               speakeriden(num_subject)=speakeriden(num_subject)+1;
            end
            a=logsig([y1 1]*w1);
            y2=logsig([a 1]*w2+[b 1]*u2);
            z=logsig([y2 1]*v2);
            z2=logsig([z 1]*v3);
            z3=[z2 1]*v4;
            Z3(p,:)=z3;
            [a1,b1]=max(z3);
            [a2,b2]=max(x0);
            if b1==b2
                subiden(label(p))=subiden(label(p))+1;
            end
            Ephonem(label(p))=sum((z3-x0).^2)+ Ephonem(label(p));
            num(label(p))=num(label(p))+1;
            Espeaker(label(p))=sum((d-b).^2) +Espeaker(label(p));
        end
    end
    speakeriden(num_subject)= (speakeriden(num_subject)./L)*100;
end
Ephonem= Ephonem./num;
Espeaker=Espeaker./num;
subiden=(subiden./num)*100;








%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eta=0.5;
max_error=0.01;
alpha=0.7;
p=1;
k=1;
E=0;
S=0;
n_z3=18;
n_z2=18;
n_z=250;
n_y2=160;
n_a=45;
n_b=100;
n_y1=160;
n_x=250;
n_x0=18;
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
u2=0.1*rand(n_b,n_y2);
du2=0.001*rand(n_b,n_y2);
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
        y2=logsig([a 1]*w2+[d]*u2);% d or b!!??
        z=logsig([y2 1]*v2);
        z2=logsig([z 1]*v3);
        z3=[z2 1]*v4;
        %%%
        E=E+sum((x0-z3).^2)+sum((d-b).^2);
        Espeach=sum((x0-z3).^2)+Espeach;
        Espeaker=Espeaker+sum((d-b).^2);
        %AA=(+0.5/500)*k+(0.01-0.5/500);
        AA=0.01;
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
        D=[d];
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

    %save('C:\Users\Admin\BME\project arshad\result\NN 8 layer-hazfe sokot\bozorgtar.mat')

    %hold on;
    %plot(e);
    %drawnow;
end




