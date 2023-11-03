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
for in=100:199
    E1=0; E2=0; PP=0;
    for out=100:199
        out;
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
                dz3=(zout-z3)*0.01;%*0.01;
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
    end
     EE=sqrt(E1/(PP*18))
        e1(k)=EE;
        EE=sqrt(E2/(PP*18))
        e2(k)=EE;
        k=k+1;
        eta=0.01;

    save('C:\Users\Admin\BME\project arshad\result\NN 8 layer-hazfe sokot\vc2')
end


