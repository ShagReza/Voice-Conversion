label=load(['C:\Users\Admin\BME\project arshad\data\LABLES\Zt',num2str(4),num2str(100)]);
label=label.Z;
j=0; vakdar400=[]; o=0; speech400=[];
i=1;
while i<=length(label)
    if label(i)<14 || label(i)==23 ||label(i)==28 || label(i)==43
        j=j+1;
        vakdar400(j)=i;
    end
    if label(i)~=40
        o=o+1;
        speech400(o)=i;
    end
    i=i+1;
end

label=load(['C:\Users\Admin\BME\project arshad\data\LABLES\Zt',num2str(5),num2str(100)]);
label=label.Z;
j=0; vakdar405=[]; i=1; speech405=[]; o=0;
while i<=length(label)
    if label(i)<14 || label(i)==23 ||label(i)==28 || label(i)==43
        j=j+1;
        vakdar405(j)=i;
    end
    if label(i)~=40
        o=o+1;
        speech405(o)=i;
    end
    i=i+1;
end


%%%
for bbb=1:1000
for out=100:199
    in=out
    for num_session=4:5
        vakdar=[]; speech=[];
        if num_session==4
            vakdar=vakdar400;
            speech =speech400;
        else
            vakdar=vakdar405;
            speech=speech405;
        end
        xx=[]; dd=[]; m=0; xx1=[]; dd1=[];
        if num_session==4 && in==168 in=in+1; end
        xx=load(['C:\Users\Admin\BME\project arshad\data\LHCB_equal ROS WAVES\',num2str(num_session),num2str(in),'_to_100']);
        xx=xx.CBw;  xx(:,19:end)=[];
        dd(1:length(xx),1:(subend-subfirst+1))=0.01;
        for i=1:length(vakdar)
            dd(vakdar(i),in-99)=.99;
        end
        j=0;
        for i=1:length(speech)
            j=j+1;
            xx1(j,:)=xx(speech(i),:);
            dd1(j,:)=dd(speech(i),:);
        end
        XX_in=xx1; DD_in=dd1;

        xx=[]; dd=[]; m=0; xx1=[]; dd1=[];
        if num_session==4 && out==168 out=out+1; end
        xx=load(['C:\Users\Admin\BME\project arshad\data\LHCB_equal ROS WAVES\',num2str(num_session),num2str(out),'_to_100']);
        xx=xx.CBw;  xx(:,19:end)=[];
        dd(1:length(xx),1:(subend-subfirst+1))=0.01;
        for i=1:length(vakdar)
            dd(vakdar(i),out-99)=0.99;
        end
        j=0;
        for i=1:length(speech)
            j=j+1;
            xx1(j,:)=xx(speech(i),:);
             dd1(j,:)=dd(speech(i),:);
        end
        XX_out=xx1; DD_out=dd1;
        %______________________________________________________________________
        P=length(XX_out);
        %r=randperm(P);
        r=1:P;
        xx_in = XX_in(r,:);  xx_out = XX_out(r,:);
        dd_in = DD_in(r,:);  dd_out = DD_out(r,:);
        for p=1:P
            x0=xx_in(p,:);
            zout=xx_out(p,:);
            d_in=dd_in(p,:);
            d_out=dd_out(p,:);
            x=logsig([x0 1]*v0);
            y1=logsig([x 1]*v1);
            b=logsig([y1 1]*u1);
            a=logsig([y1 1]*w1);
            y2=logsig([a 1]*w2+[d_out 1]*u2);
            z=logsig([y2 1]*v2);
            z2=logsig([z 1]*v3);
            z3=[z2 1]*v4;
            %%
            E1=E1+sum((zout-z3).^2);
            E2=E2+sum((d_in-b).^2);
            dz3=(zout-z3)*1.2;%*0.01;
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
            D=[d_out 1];
            du2=eta.*D'*dY2+alpha.*du2;
            u2=u2+du2;
            db=(d_in-b).*b.*(1-b);
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
        PP=P+PP;
    end
    EE=sqrt(E1/(PP*18))
    e1(k)=EE;
    EE=sqrt(E2/(PP*18))
    e2(k)=EE;
    k=k+1;
    eta=0.05;

end
save('C:\Users\Admin\BME\project arshad\result\NN 8 layer-hazfe sokot\vc3')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% tabdile seda ba shabakeye khodanjomani %%%
num_session=4
marja=100;
hadaf=102;
num_subject=marja;
xx=[]; dd=[];
xx=load(['C:\Users\Admin\BME\project arshad\data\LHCB_equal ROS WAVES\',num2str(num_session),num2str(marja),'_to_100']);
xx=xx.CBw;  xx(:,19:end)=[];

num_subject=hadaf;
xxhadaf=[]; dd=[];
xxhadaf=load(['C:\Users\Admin\BME\project arshad\data\LHCB_equal ROS WAVES\',num2str(num_session),num2str(hadaf),'_to_100']);
xxhadaf=xxhadaf.CBw; xxhadaf(:,19:175)=[];
dd(1:length(xxhadaf),1:(subend-subfirst+1))=0.01;
dd(1:length(xxhadaf),num_subject-99)=.99; %%cod goyande az 1 be .95

B=[]; A=[]; Z3=[]; B=[];
P=length(xx);
for p=1:P
    x0=xx(p,:);
    d=dd(p,:);
    x=logsig([x0 1]*v0);
    y1=logsig([x 1]*v1);
    b=logsig([y1 1]*u1);
    B(p,:)=b;
    a=logsig([y1 1]*w1);
    y2=logsig([a 1]*w2+[d 1]*u2);
    z=logsig([y2 1]*v2);
    z2=logsig([z 1]*v3);
    z3=[z2 1]*v4;
    Z3(p,:)=z3;
end
%plot(xx(:,2)),hold on,plot(Z3(:,2),'r')

% [a,b]=max(B');
% plot(b,'.');
%%%%%%%%%%%%%%%% SNR %%%%%%%%%%%%%%%
% noise=Z3-xx;
% n=noise.^2;
% pownoise=sum(sum(n));
% powsignal=sum(sum(xx.^2));
% SNR=powsignal/pownoise;
% l=l+1;
% snrtest(l)=SNR;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CBw=[];
CBw(:,1:18)=Z3(:,1:18);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_CB=load('C:\Users\Admin\BME\project arshad\programs\p5_mean & var extraction\mean_CB.mat');
mean_CB=mean_CB.mean_CB;
var_CB=load('C:\Users\Admin\BME\project arshad\programs\p5_mean & var extraction\var_CB.mat');
var_CB=var_CB.var_CB;
num_subject=marja;
mean_CB1=mean_CB(num_subject,:);
var_CB1=var_CB(num_subject,:);

CBw(:,1:18) = CBw(:,1:18) .* repmat(var_CB1, size(CBw , 1) ,1) ;
CBw(:,1:18) = CBw(:,1:18) + repmat(mean_CB1 , size(CBw , 1),1) ;
%%%%%%%%%% bazsazi soot az khoroji shabake %%%%%%%%%%
%pitch:
filename = ['C:\Users\Admin\BME\project arshad\data\equal ROS waves\W',num2str(num_session),num2str(num_subject),'_to_100.wav'];
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
if Pitch(2)==0
    Pitch(2)=1;
end

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

wavplay(smp*10,fs)
%filename=['C:\Users\Admin\BME\project arshad\gozaresh prozhe\8layer\',num2str(ns),'prime.wav'];
%wavwrite(smp,fs,filename);


