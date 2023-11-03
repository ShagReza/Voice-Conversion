%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% moshref be hadaf (nural net) %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;  clear all; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load
input=FFT175;
size_input = size(input);
P=size_input(1);
%%%
eta=.1;
MaxE=.01;
alpha=.7;
M=175;   %first layer
N=256;    %second layer
%%%
Wa=0.01*rand(M+1,N)-0.005;
Vb=0.01*rand(N+1,M)-0.005;
V=0.01*rand(M,M)-0.005;

DeltaWa=zeros(M+1,N);
DeltaVb=zeros(N+1,M);
DeltaV=zeros(M,M);


clear CB r size_input
p=1;k=1;E=0;counter=1;F=1;q=1;
tic
%%%%
while (F)
    r=randperm(P);
    input = input(r,:);
    d = input;%desired output
    while p<=P
        xa=input(p,:);
        xa=[xa  1];
        ya=logsig(xa*Wa);
        ya = [ya , +1];
        xb=logsig(ya*Vb);
        x=(xb*V);
        %E=E+sum((d(p,:)-xb).^2);
        E=E+sum((d(p,:)-x).^2);
        
        
        dx=(d(p,:)-x)*.001;
        dxb=(dx*V') .* xb .* (1-xb);
        %dxb=(d(p,:)-xb).* xb .* (1-xb);
        dya=(dxb*Vb') .* ya .* (1-ya);
        dya(end) = [];
        DeltaV=eta.* xb' *dx+alpha.*DeltaV;
        V=V+DeltaV;
        DeltaVb=eta.* ya' *dxb+alpha.*DeltaVb;
        Vb=Vb+DeltaVb;
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
end
