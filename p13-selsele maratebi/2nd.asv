load
Wa1=Wa; Va1=Va;

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
%eta=(-0.4/4999)*k+(0.5+0.4/4999);
while (F) && (eta>0)
    k
    r=randperm(P);
    input = input(r,:);
    d = input;   %desired output
    while p<=P
        xa=input(p,:);
        xa=[xa  1];
        ya1=logsig(xa*Wa1);
        ya1 = [ya1 , +1];
        z1=logsig(ya1*Va1);
        ya=logsig(xa*Wa);
        ya = [ya , +1];
        z=logsig(ya*Va);
        z = [z , z1];
        yb=logsig(z*Wb);
        yb = [yb , +1];
        xb=(yb*Vb);

        E=E+sum((d(p,:)-xb).^2);

        dxb=(d(p,:)-xb)*.01;
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
    %eta=(-0.4/4999)*k+(0.5+0.4/4999);
    eta=0.01;
end


% ------loading data for test----------------------
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
    ya1=logsig(xa*Wa1);
    ya1 = [ya1 , +1];
    z1=logsig(ya1*Va1);
    ya=logsig(xa*Wa);
    ya = [ya , +1];
    z=logsig(ya*Va);
    z = [z , z1];
    yb=logsig(z*Wb);
    yb = [yb , +1];
    xb=(yb*Vb);
    XB(p,:)=xb;
    p=p+1;
end





input=xx;
size_select = size(input);
P=size_select(1);
p=1;
XB=[]; Z=[];
while p<=P
    xa=input(p,:);
    xa=[xa  1];
    ya1=logsig(xa*Wa1);
    ya1 = [ya1 , +1];
    z1=logsig(ya1*Va1);
    Z1(p)=z1;
    ya=logsig(xa*Wa);
    ya = [ya , +1];
    z=logsig(ya*Va);
    z = [z z1];
    %Z(p)=z;
    yb=logsig(z*Wb);
    yb = [yb , +1];
    xb=(yb*Vb);
    XB(p,:)=xb;
    p=p+1;
end


WBB=Wb;
input=xx;
size_select = size(input);
P=size_select(1);
p=1;
XB=[]; Z=[];
Wb(2,:)=[];
while p<=P
    xa=input(p,:);
    xa=[xa  1];
    ya1=logsig(xa*Wa1);
    ya1 = [ya1 , +1];
    z1=logsig(ya1*Va1);
    Z1(p)=z1;
    ya=logsig(xa*Wa);
    ya = [ya , +1];
    z=logsig(ya*Va);
    %z = [z z1];
    %Z(p)=z;
    
    yb=logsig(z*Wb);
    yb = [yb , +1];
    xb=(yb*Vb);
    XB(p,:)=xb;
    p=p+1;
end
subplot(1,2,1),plot(xx(:,3),xx(:,4),'.r'),hold on,plot(XB(:,3),XB(:,4),'.b'),