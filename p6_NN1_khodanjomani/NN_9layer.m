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
%takhsise code 0 be sokot:
P=length(dd);
for i=1:P
    if totallabel(i)==40
        [maxdd,indmax_dd]=max(dd(i,:));
        dd(i,indmax_dd)=.01;
    end
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
n_z3=18;
n_z2=60;
n_z=100;
n_y2=60;
n_a=45;
n_b=30;
n_y1=60;
n_x=100;
n_x2=60;
n_x0=18;
%%%
v0=0.1*rand(n_x0+1,n_x2);
dv0=0.001*rand(n_x0+1,n_x2);
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
v5=0.1*rand(n_x2+1,n_x);
dv5=0.001*rand(n_x2+1,n_x);
E=5; k=1; e=[];
%%%
P=length(xx);
eta=(-0.4/500)*k+(0.5+0.4/500);
while E>0.01
    E=0;
    r=randperm(P);
    xx = xx(r,:);
    dd = dd(r,:);
    for p=1:P
        x0=xx(p,:);
        d=dd(p,:);
        x2=logsig([x0 1]*v0);
        x=logsig([x2 1]*v5);
        y1=logsig([x 1]*v1);
        b=logsig([y1 1]*u1);
        a=logsig([y1 1]*w1);
        y2=logsig([a 1]*w2+[b 1]*u2);
        z=logsig([y2 1]*v2);
        z2=logsig([z 1]*v3);
        z3=[z2 1]*v4;
        %%%
        E=E+sum((x0-z3).^2)+sum((d-b).^2);
        dz3=(x0-z3)*0.01;
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
        dY1(n_y1+1)=[];
        X=[x 1];
        dv1=eta.*X'*dY1+alpha.*dv1;
        v1=v1+dv1;
        dX=X.*(1-X).*(dY1*v1');
        dX(n_x+1)=[];
        X2=[x2 1];
        dv5=eta.*X2'*dX+alpha.*dv5;
        v5=v5+dv5;
        dX2=X2.*(1-X2).*(dX*v5');
        dX2(n_x2+1)=[];
        X0=[x0 1];
        dv0=eta.*X0'*dX2+alpha.*dv0;
        v0=v0+dv0;
    end
    EE=sqrt(E/(length(xx)*18))
    e(k)=EE;
    k=k+1
    if k<600
        eta=(-0.4/500)*k+(0.5+0.4/500)
    else
        eta=0.01;
    end
    if mod(k,20)==0
        %save('C:\Users\Admin\BME\project arshad\result\NN_8layer(18-100-60-(5-30)-60-100-18)')
    end
    %hold on;
    %plot(e);
    %drawnow;
end
%%%%%%%%%%%%%%%%%%%%