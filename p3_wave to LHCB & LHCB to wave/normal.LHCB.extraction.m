%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% wave to LHCB_normalized &  fs=22050 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; close all; clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ROS signal dovom manand signal aval mishaad
mean_CB=load('C:\Users\Admin\BME\project arshad\programs\p5_mean & var extraction\mean_CB.mat');
mean_CB=mean_CB.mean_CB;
var_CB=load('C:\Users\Admin\BME\project arshad\programs\p5_mean & var extraction\var_CB.mat');
var_CB=var_CB.var_CB;

for num_session=4:5
    for num_subject=65:303 
        num_subject
        load( ['C:\Users\Admin\BME\project arshad\data\LHCB_without normalization\LHCB_',num2str(num_session),num2str(num_subject)],'CBw' );
        %hazfe miyangin o variyans az LHCB:
        mean_CBw=mean_CB(num_subject,:);
        var_CBw=var_CB(num_subject,:);
        CBw(:,1:18) = CBw(:,1:18) - repmat(mean_CBw , size(CBw , 1) ,1) ;
        CBw(:,1:18) = CBw(:,1:18) ./ repmat(var_CBw , size(CBw , 1) ,1) ;
        save( ['C:\Users\Admin\BME\project arshad\data\LHCB_ normalized\LHCB_',num2str(num_session),num2str(num_subject)],'CBw' );
    end
end
