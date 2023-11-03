clear all,close all,clc

for num_subject=100:130 %wav W564 nadarim
        num_subject
        CB=load(['C:\Users\Admin\BME\project arshad\data\LHCB_without normalization\LHCB_',num2str(5),num2str(100)]);
        CB=CB.CBw;
        filename = ['C:\Users\Admin\BME\project arshad\data\equal ROS waves\W',num2str(5),num2str(num_subject),'_to_100.wav'];
        [s1,fs] = wavread(filename);
        %%%
        %-feature extraction--------------------------------------------------
        N = 512 ; %1024		            % Frame length
        stp = 256;	%512	            % Frame rate
        L = 18 ;		            % Number of LHCB parameters
        Q = 18 ;		            % Number of Filters
        fc = 7500 ;	                % Final frequency is 7500 Hz
        fs = 22050 ;  %44100          	% Sampling frequency
        kc = round( fc * N / fs );	% Equivallent of 7500 Hz in DFT domain
        f = fs / N * (0:kc);		% Frequencies samples in hertz
        win =[0.54+0.46*cos(2*pi/511*(-256:255))]';	   % Temporal window
        fbark = 6 * log( f / 600 + sqrt( ( f / 600 ).^2 + 1 ) ) ;
        CBFBS1 = zeros(Q,length(f));

        for zk = 1 : Q
            bk = ( ( fbark >= zk - 1 ) & ( fbark < zk +1 ) ) ;
            CBFBS1(zk , :) = bk .* ( 0.5 + 0.5 * cos( pi * ( fbark - zk ) ) ) ;
        end

        CBFBS1 = CBFBS1'; CBFBS2 = CBFBS1 .^ 2 ;
        CBFBS2(1:175,19:175)=rand(175,157); CBFBS=CBFBS2;
        CBw = []; CB2 = [];
        m = 1; n=1; x=1;
        cbfinal1=[]; cbfinal2=[];
        while m < length(s1) - N + 1,
            frm = s1 ( m : m + N - 1 );
            frm = frm - mean(frm);
            frmwin = frm .* win;
            frmfft = fft ( frmwin , N ) ;
            frmfft ( kc+2 : end ) = [];
            frmfftsqr =(( real ( frmfft .* conj ( frmfft ) ))*4)';
            %energy of filter
            cb = frmfftsqr * CBFBS ;
            cbfinal1=[cbfinal1 cb];
            cbcepst = log ( 0.01 + cb ) ;
            cbfinal2=[cbfinal2 cbcepst];
            %after log
            CBw = [CBw ; cbcepst] ;
            m = m + 256 ;  n = n + 1;
        end
        CB=CB(1:459,:);
        CBw=CBw(1:459,:);
        score_GA(num_subject)=sum(sum((CBw-CB).^2))/length(CBw);
        clear CB; clear CBw
end
    

%-----
for num_subject=100:120 %wav W564 nadarim
        CB=load(['C:\Users\Admin\BME\project arshad\data\LHCB_without normalization\LHCB_',num2str(5),num2str(100)]);
        CB=CB.CBw;
        filename = ['C:\Users\Admin\BME\project arshad\data\equal ROS waves\W',num2str(5),num2str(num_subject),'_DTW_100.wav'];
        [s1,fs] = wavread(filename);
        %%%
        %-feature extraction--------------------------------------------------
        N = 512 ; %1024		            % Frame length
        stp = 256;	%512	            % Frame rate
        L = 18 ;		            % Number of LHCB parameters
        Q = 18 ;		            % Number of Filters
        fc = 7500 ;	                % Final frequency is 7500 Hz
        fs = 22050 ;  %44100          	% Sampling frequency
        kc = round( fc * N / fs );	% Equivallent of 7500 Hz in DFT domain
        f = fs / N * (0:kc);		% Frequencies samples in hertz
        win =[0.54+0.46*cos(2*pi/511*(-256:255))]';	   % Temporal window
        fbark = 6 * log( f / 600 + sqrt( ( f / 600 ).^2 + 1 ) ) ;
        CBFBS1 = zeros(Q,length(f));

        for zk = 1 : Q
            bk = ( ( fbark >= zk - 1 ) & ( fbark < zk +1 ) ) ;
            CBFBS1(zk , :) = bk .* ( 0.5 + 0.5 * cos( pi * ( fbark - zk ) ) ) ;
        end

        CBFBS1 = CBFBS1'; CBFBS2 = CBFBS1 .^ 2 ;
        CBFBS2(1:175,19:175)=rand(175,157); CBFBS=CBFBS2;
        CBw = []; CB2 = [];
        m = 1; n=1; x=1;
        cbfinal1=[]; cbfinal2=[];
        while m < length(s1) - N + 1,
            frm = s1 ( m : m + N - 1 );
            frm = frm - mean(frm);
            frmwin = frm .* win;
            frmfft = fft ( frmwin , N ) ;
            frmfft ( kc+2 : end ) = [];
            frmfftsqr =(( real ( frmfft .* conj ( frmfft ) ))*4)';
            %energy of filter
            cb = frmfftsqr * CBFBS ;
            cbfinal1=[cbfinal1 cb];
            cbcepst = log ( 0.01 + cb ) ;
            cbfinal2=[cbfinal2 cbcepst];
            %after log
            CBw = [CBw ; cbcepst] ;
            m = m + 256 ;  n = n + 1;
        end
        CB=CB(1:459,:);
        CBw=CBw(1:459,:);
        score_DTW(num_subject)=sum(sum((CBw-CB).^2))/length(CBw);
        clear CB;
end
  plot(t,score_DTW(100:105)),hold on,plot(t,score_GA(100:105),'r')
  %%%%%

 plot(CB(n,1:18)),hold on,plot(CBw(n,1:18),'r'),hold on,plot(CBw1(n,1:18),'g')
 
 f=[]; g=[];
 for i=1:459
     f(i)=sum((CBw1(i,1:18)-CB(i,1:18)).^2);
     g(i)=sum((CBw(i,1:18)-CB(i,1:18)).^2);
 end
     
 
 
 %%%%
size1=size(CBw1);
size2=size(CBw);
size3=size(CB);
delta2=[];
delta1=[];
delta3=[];
for i=1:(size1(1)-1)
    delta1=[delta1,abs(mean(CBw1(i,:)-CBw1(i+1,:)))];
end
for i=1:(size2(1)-1)
    delta2=[delta2,abs(mean(CBw(i,:)-CBw(i+1,:)))];
end
for i=1:(size3(1)-1)
    delta3=[delta3,abs(mean(CB(i,:)-CB(i+1,:)))];
end
plot(delta3),hold on,plot(delta1,'g'),hold on,plot(delta2,'r')

