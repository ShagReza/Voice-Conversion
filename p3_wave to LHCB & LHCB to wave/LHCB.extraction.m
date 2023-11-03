%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% wave to LHCB_without normalization &  fs=22050 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;  clear all; close all

for num_session=4
    for num_subject=100:199 %wav W564 nadarim
        num_subject
        filename = ['C:\Users\Admin\BME\project arshad\data\WAV\W',num2str(num_session),num2str(num_subject),'.wav'];
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
        save( ['C:\Users\Admin\BME\project arshad\data\LHCB_without normalization\LHCB_',num2str(num_session),num2str(num_subject)],'CBw' );
        clear CB;
    end
end