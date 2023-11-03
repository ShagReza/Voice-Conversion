%load
v1=0.01*rand(128+1,50)-0.005;
v2=0.01*rand(50+1,5)-0.005;
v3=0.01*rand(5,50)-0.005;
v4=0.01*rand(50+1,128)-0.005;

Deltav1=zeros(128+1,50);
Deltav2=zeros(50+1,5);
Deltav3=zeros(5,50);
Deltav4=zeros(50+1,128);


while (F) && (eta>0)
    k
    r=randperm(P);
    input = input(r,:);
    d = input;   %desired output
    while p<=P
        xa=input(p,:);
        D=dd(p,:);
        xa=[xa  1];
        ya=logsig(xa*Wa);
        ya = [ya , +1];
        z=logsig(ya*Va);
        b1=logsig(ya*v1);
        b1=[b1 , 1];
        b2=logsig(b1*v2);
        b3=logsig(D*v3);
        b3=[b3 ,1];
        yb=logsig(b3*v4 + z*Wb);
        yb = [yb , +1];
        xb=(yb*Vb);
       
        E=E+sum((d(p,:)-xb).^2);

        dxb=(d(p,:)-xb)*.001;
        dyb=(dxb*Vb') .* yb .* (1-yb);
        dyb(end) = [];
        db3=(dyb*v4') .* b3 .* (1-b3);
        db3(end) = [];
        db2=(D-b2).*b2.*(1-b2);
        db1=(db2*v2') .* b1 .* (1-b1);
        db1(end) = [];
       

        Deltav4=eta.* b3' *dyb+alpha.*Deltav4;
        v4=v4+Deltav4;
        Deltav3=eta.* b2' *db3+alpha.*Deltav3;
        v3=v3+Deltav3;
        Deltav2=eta.* b1' *db2+alpha.*Deltav2;
        v2=v2+Deltav2;
        Deltav1=eta.* ya' *db1+alpha.*Deltav1;
        v1=v1+Deltav1;
        
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