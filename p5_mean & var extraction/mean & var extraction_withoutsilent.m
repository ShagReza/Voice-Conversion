%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% In the name of God %%%%%%%
%%%%%%% shaghayegh reza %%%%%%%%%
%%%%%%%% project arshad %%%%%%%%%
%%%mean & var extraction %%%%%%%%
%%%%%%%%%% 25 ordibehesht 87 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; close all; clear all;
%dar in barname dadegane 100 nafar do jomle bedone sokot normalise mishavad
%meghdar mean va var morede estefade baraye normalize baraye hame yeksan
%ast va dar matris dd ham kode goyande gharar migirad
xx=[];
dd=[];
m=0;
num_session=4;
subfirst=100;
subend=199;
for num_subject=subfirst:subend
    for num_session=4:5
        x=load(['C:\Users\Admin\BME\project arshad\data\LHCB_without normalization\LHCB_',num2str(num_session),num2str(num_subject)]);
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
%hazfe sokot:
P=length(dd);
i=1;
while i < P
    if totallabel(i)==40
        xx(i,:)=[];
        dd(i,:)=[];
        totallabel(i)=[];
        i=i-1;
        P=P-1;
    end
    i=i+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CB=xx;
mean_CB=mean(CB);
var_CB=sqrt(sum((CB-repmat(mean_CB , size(CB , 1) , 1)).^2)/length(CB));
%%%%%normalize kardane parametrhaye LHCB:
CB(:,1:18) = CB(:,1:18) - repmat(mean_CB , size(CB , 1) ,1) ;
CB(:,1:18) = CB(:,1:18) ./ repmat(var_CB , size(CB , 1) ,1) ;

XX=xx;
xx=CB;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%normalize be sokot va baraye har fard mojaza:

clc; close all; clear all;
%dar in barname dadegane 100 nafar do jomle bedone sokot normalise mishavad
%meghdar mean va var morede estefade baraye normalize baraye hame yeksan
%ast va dar matris dd ham kode goyande gharar migirad
xx=[];
dd=[];
m=0;
num_session=4;
subfirst=100;
subend=199;
totallabel=[];
Label=[];
XX=[];
DD=[];

for num_subject=subfirst:subend
    xx=[]; x=[];
    Label=[]; label=[];
    dd=[]; m=0; d=[];
    for num_session=4:5
        x=load(['C:\Users\Admin\BME\project arshad\data\LHCB_without normalization\LHCB_',num2str(num_session),num2str(num_subject)]);
        x=x.CBw;
        xx=[xx;x];
        d(1+m:length(x)+m,1:(subend-subfirst+1))=0.01;
        d(1+m:length(x)+m,num_subject-99)=.99; %%cod goyande az 1 be .95
        dd=d;
        m=m+length(x);
        label=load(['C:\Users\Admin\BME\project arshad\data\LABLES\Zt',num2str(num_session),num2str(num_subject)]);
        label=label.Z;
        Label=[Label;label];
        totallabel=[totallabel;label];
    end
    
    P=length(xx);
    i=1;
    while i < P+1
        if Label(i)==40
            xx(i,:)=[];
            dd(i,:)=[];
            Label(i)=[];
            i=i-1;
            P=P-1;
        end
        i=i+1;
    end
    DD=[DD;dd];
    xx(:,19:end)=[];
    mean_xx(num_subject,1:18)=mean(xx);
    var_xx(num_subject,1:18)=sqrt(sum((xx-repmat(mean_xx(num_subject,1:18) , size(xx , 1) , 1)).^2)/length(xx));

    xx(:,1:18) = xx(:,1:18) - repmat(mean_xx(num_subject,1:18) , size(xx , 1) ,1) ;
    xx(:,1:18) = xx(:,1:18) ./ repmat(var_xx(num_subject,1:18) , size(xx , 1) ,1) ;
    XX=[XX;xx];
end
xx=[]; xx=XX;
dd=[]; dd=DD; 
