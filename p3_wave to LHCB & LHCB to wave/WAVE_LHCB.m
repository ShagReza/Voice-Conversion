%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% wave to LHCB & LHCB to wave %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;  clear all; close all
N=1024; n=5; K=105;%184
[s1,fs] = wavread ( 'C:\Users\Admin\BME\project arshad\data\Wav\W4100.WAV');

%%%%%%%%%%%%%%%%%%%%%%%%% WAVE to LHCB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------
%   Step 1 :  Pitch Extraction
%------------------------------------------------

m=1;  Pitch=[];  PP=0;
while m < length(s1) - N + 1,
    frm = s1( m : m + N - 1 );
    a = real ( ifft ( abs ( fft ( frm , N ) ) ) );
    [cc,bb] = max ( a ( 100 : 512 ) ) ;
    pitch = bb + 99 ;
    PP = PP + pitch ;
    m = m + pitch;
    Pitch=[Pitch pitch];
end
m=1; i=1;
%%%
while i <= length ( Pitch )
    m = m + Pitch ( i )  ;
    if m >= 512
        for j = length ( Pitch ) + 1 : -1 : i +1
            Pitch ( j ) = Pitch ( j - 1 ) ;
        end
        Pitch ( i ) = 512;
        Pitch ( i + 1 ) = m - 512 ;
        m = 0;
    end
    i = i + 1;
end

%-----------------------------------------------
%   Step 2 :  LHCB Extraction
%------------------------------------------------

N = 1024 ;							% Frame length
stp = 512;							% Frame rate
L = 18 ;							% Number of LHCB parameters
Q = 18 ;							% Number of Filters
fc = 7500 ;							% Final frequency is 7500 Hz
fs = 44100 ;						% Sampling frequency
kc = round( fc * N / fs );			% Equivallent of 7500 Hz in DFT domain
f = fs / N * (0:kc);				% Frequencies samples in hertz
win = hamming(N);					% Temporal window
fbark = 6 * log( f / 600 + sqrt( ( f / 600 ).^2 + 1 ) ) ;
CBFBS = zeros(Q,length(f));

for zk = 1 : Q
   bk = ( ( fbark >= zk - 1 ) & ( fbark < zk +1 ) ) ;
   CBFBS(zk , :) = bk .* ( 0.5 + 0.5 * cos( pi * ( fbark - zk ) ) ) ;
end

CBFBS = CBFBS';
CBFBS = CBFBS .^ 2 ;
CBFBS(1:175,19:175)=rand(175,157);

CBw = []; CB2 = [];
m = 1; n=1; x=1;
cbfinal1=[]; cbfinal2=[];
while m < length(s1) - N + 1,
    frm = s1 ( m : m + N - 1 );
    frm = frm - mean(frm);
    frmwin = frm .* win;
    frmfft = fft ( frmwin , N ) ;
    frmfft ( kc+2 : end ) = [];
    frmfftsqr = real ( frmfft .* conj ( frmfft ) )';
    cb = frmfftsqr * CBFBS ;%energy of filter
    cbcepst = log ( 0.01 + cb ) ;
    CBw = [CBw ; cbcepst] ;%after log
    m = m + 512 ;
    n = n + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% LHCB to WAVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m1=1; N2=175; m2=1;
newfrms = exp ( CBw ) - 0.01 ;
z=size(newfrms,1)*1024-(size(newfrms,1)-1)*512;
smp=zeros(1,z);
PP=0;  q=0;
  
kk = newfrms( m2,:) / CBFBS;
kk ( 1024 : -1 : 851 ) = kk ( 2 : 175 );  % mirrorong
kk ( 176 : 850 ) = 0;
dd =  real ( ifft ( sqrt(kk) , 1024 ) ) ;
d = dd .* hamming(1024)' ;  %for deletting the effect of needle picks at first and end of frames
d2 = [d(513:1024),d(1:512)];
smp ( 512 * q + 1 : 512 * q + 1024 ) = smp ( 512 * q + 1 : 512 * q + 1024 ) + d2 ;

while ( Pitch ( m1 ) ~= 512 ) & ( m2 == 1 )
    m1=m1+1;
end

m2 = m2 + 1; m1 = m1 + 1;

while ( m2 <= size(newfrms,1) ) & ( m1 < size ( Pitch , 2 ) )

    kk = newfrms( m2,:) / CBFBS;
    kk ( 1024 : -1 : 851 ) = kk ( 2 : 175 );  % mirroring
    kk ( 176 : 850 ) = 0;
    dd =  real ( ifft ( sqrt(kk) , 1024 ) ) ;
    d = dd .* hamming ( 1024 )' ;  %for deletting the effect of needle picks at first and end of frames
    d2 = [d(513:1024),d(1:512)];

    PP = Pitch ( m1 );

    while ( Pitch ( m1 ) ~= 512 ) & ( m1 < size ( Pitch , 2 ) )
        smp ( 512 * q + PP : 512 * q + 1024 + PP  -1) = smp ( 512 * q + PP : 512 * q + 1024 + PP  -1) + d2 ;
        m1 = m1 + 1 ;
        PP = PP + Pitch ( m1 ) -1 ;
    end

    q = q + 1;
    m2 = m2 + 1;
    m1 = m1 + 1;
    newfrms11=newfrms(:,1:18).^2;
    newfrms1=sum( newfrms11);
    d211=d2.^2;
    d21=sum(d211);
end

wavwrite ( smp, fs ,'C:\Users\Admin\BME\project arshad\data\Wav\W4100_2.WAV' ) ;
wavplay(s1,44100)
pause(1);
wavplay(smp,44100)




