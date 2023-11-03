%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% In the name of God %%%%%%%
%%%%%%% shaghayegh reza %%%%%%%%%
%%%%%%%% project arshad %%%%%%%%%
%%%mean & var extraction %%%%%%%%
%%%%%%%%%% 25 ordibehesht 87 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; close all; clear all;
%dar in barname maghadir ra baraye ferquency 22050 tanzim kardeiim:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 512;%1024							 % Frame length
stp = 256;%512							 % Frame rate
L = 18 ;								 % Number of LHCB parameters
Q = 18 ;								 % Number of Filters
fc = 7500 ;								 % Final frequency is 7500 Hz
fs = 22050;%44100							 % Sampling frequency
kc = round( fc * N / fs );				 % Equivallent of 7500 Hz in DFT domain
f = fs / N * (0:kc);					 % Frequencies samples in hertz
win=[0.54+0.46*cos(2*pi/511*(-256:255))]'; % Temporal window
fbark = 6 * log( f / 600 + sqrt( ( f / 600 ).^2 + 1 ) ) ;
CBFBS = zeros(Q,length(f));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for zk = 1 : Q
    bk = ( ( fbark >= zk - 1 ) & ( fbark < zk +1 ) ) ;
    CBFBS(zk , :) = bk .* ( 0.5 + 0.5 * cos( pi * ( fbark - zk ) ) ) ;
end
CBFBS = CBFBS';
CBFBS = CBFBS .^ 2 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for K=1:304
    K
    filename = ['C:\Users\Admin\BME\project arshad\data\session2_farsdot\S2',num2str(K),'.wav'];
    [s,fs] = wavread(filename);
    CB = []; m = 1;
    while m < length(s) - N + 1,
        frm = s(m : m+N-1);
        frm = frm - mean(frm);
        frmwin = frm .* win;
        frmfft = fft(frmwin , N);
        frmfft(kc+2 : end) = [];
        frmfftsqr = [real( frmfft .* conj(frmfft) )*4]';
        cb = frmfftsqr * CBFBS ;
        cbcepst = log(0.01 + cb) ;
        CB = [CB ; cbcepst] ;
        m = m + stp;
    end
    mean_CB(K,:)=mean(CB);
    var_CB(K,:)=sqrt(sum((CB-repmat(mean_CB(K,:) , size(CB , 1) , 1)).^2)/length(CB));
    %var_CB(K,:)=sqrt(var(CB,1)); %eshtebah ast
    clear CB;
end

