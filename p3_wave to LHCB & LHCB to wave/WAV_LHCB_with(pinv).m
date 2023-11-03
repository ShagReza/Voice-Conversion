%  In the name of god
%84/2/25
% This program extracts LHCB parameters with gaame pishravi equal to
% 512 for all of the frames, but speech reconstructio is by gaame pishravi
% equal to Pitch for each frame.
clc
clear all
close all

%   Step 1 : LHCB Extraction
%   ========================

%    In this part parameters are extracted with gaame pishravi=512

filename = ['C:\Users\Admin\BME\project arshad\data\Wave\W4101.WAV'];
[s,fs] = wavread ( filename );

N = 1024 ;										% Frame length
stp = N/2;										% Frame rate

L = 18 ;											% Number of LHCB parameters
Q = 18 ;											% Number of Filters
fc = 7500 ;										% Final frequency is 7500 Hz
kc = round( fc * N / fs );					% Equivallent of 7500 Hz in DFT domain  
f = fs / N * (0:kc);							% Frequencies samples in hertz

%f = fs / N * ( 0 : ( ( ( exp ( 30 ) ) ^ 2 - 1 ) * 300 ) / exp ( 30 ) );

fbark = 6 * log( f / 600 + sqrt( ( f / 600 ).^2 + 1 ) ) ;
win = hamming(N);								% Temporal window
CBFBS = zeros(Q,length(f));

for zk = 1 : Q
    bk = ( ( fbark >= zk - 1 ) & ( fbark < zk +1 ) ) ;
    CBFBS(zk , :) = bk .* ( 0.5 + 0.5 * cos( pi * ( fbark - zk ) ) ) ;
end

CBFBS = CBFBS';
CBFBS = CBFBS .^ 2 ;


% CBFBS(1:kc+1,L+1:kc+1)=0.01*rand(kc+1,kc+1-L);
% CBFBS(1:175,19:175)=rand(175,157);

CB = [];
CB2 = [];
m = 1;
n=1;

while m < length(s) - N + 1,
    frm = s ( m : m + N - 1 );
    frm = frm - mean(frm);
    frmwin = frm .* win;
    frmfft = fft ( frmwin , N ) ;
    frmfft ( kc+2 : end ) = [];
    frmfftsqr = real ( frmfft .* conj ( frmfft ) )';
    cb = frmfftsqr * CBFBS ;
    cbcepst = log ( 0.01 + cb ) ;
    CB = [CB ; cbcepst] ;
    m = m + stp ;%- 1; %Pitch ( n );
    n = n + 1;
end




%   Step 2 :  Pitch Extraction
%   ==========================


m = 1;
Pitch=[];
PP=0;
Per=100;        % This number is suitable for fs (sampling frequency) =1024

while m < length(s) - N + 1,
    
    frm = s( m : m + N - 1 );
    a = real ( ifft ( abs ( fft ( frm , N ) ) ) );
    [cc,bb] = max ( a ( Per : N/2 ) ) ;
    pitch = bb + Per-1 ;
    PP = PP + pitch ;
    m = m + pitch;
    Pitch=[Pitch pitch];
end

m=1;
i=1;

while i <= length ( Pitch ) 
    
    m = m + Pitch ( i )  ;
    if m >= N/2
        for j = length ( Pitch ) + 1 : -1 : i +1
            Pitch ( j ) = Pitch ( j - 1 ) ;
        end
        Pitch ( i ) = N/2;
        Pitch ( i + 1 ) = m - N/2 ;
        m = 0;
    end
    i = i + 1;
    
end





%   Step 3 : Normalization
%   ====================== 

%***********************
% generating K1 and K2:
% ======================

% LhcbFile = ' D:\LHCB\W' ;
% NormalizedLhcbFile = 'D:\2ndNormlhcb\' ;


% Calculation of General K1 & K2
kc=17;

K1 = zeros(1 , kc+1) ;
FileSize=0;
% for n =4:5
%     for K=100:200
%         eval(['load ',LhcbFile,num2str(n),num2str(K),'.mat' ] ) ;
%         CB2=CB(:,1:18);
        FileSize = FileSize+size(CB , 1) ;
        K1 = K1 + sum(CB , 1) ;
%     end
% end
K1 = K1 / FileSize;

K2 = zeros(1 , kc+1) ;
% for n =4:5
%     for K=100:200
%         eval(['load ',LhcbFile,num2str(n),num2str(K),'.mat' ] ) ;
%         CB2=CB(:,1:18);
        K2 = K2 + sum( ( CB - repmat(K1 , size(CB , 1) , 1) ).^2 , 1 ) ;
%     end
% end
K2 = sqrt(K2/ FileSize)  ;

%*******************************


CB(:,1:kc+1) = CB(:,1:kc+1) - repmat(K1 , size(CB(:,1:kc+1) , 1) ,1) ;
CB(:,1:kc+1) = CB(:,1:kc+1) ./ repmat(K2 , size(CB(:,1:kc+1) , 1) ,1) ;


%   Step 4 : Inverse Normalization steps 
%   ====================================

CB(:,1:kc+1) = CB(:,1:kc+1) .* repmat(K2 , size(CB(:,1:kc+1) , 1) ,1) ;
CB(:,1:kc+1) = CB(:,1:kc+1) + repmat(K1 , size(CB(:,1:kc+1) , 1) ,1) ;




%   Step 5 : Inverse LHCB steps 
%   ===========================


m1=1;
N2=kc+1;
m2=1;
newfrms = exp ( CB ) - 0.01 ;
z=size(newfrms,1)*N-(size(newfrms,1)-1)*N/2;
smp=zeros(1,z);
PP=0;
q=0;


kk = newfrms( m2,:)*pinv(CBFBS);  %/(CBFBS);
kk ( N : -1 : N-kc+1 ) = kk ( 2 : kc+1 );  %mirroring
kk ( kc+2 : N-kc ) = 0;
dd =  real ( ifft ( sqrt(kk) , N ) ) ;
d = dd .* hamming(N)' ;  %for deletting the effect of needle picks at first and end of frames 
d2 = [d(N/2+1:N),d(1:N/2)]; 
%d2=d;
smp ( N/2 * q + 1 : N/2 * q + N ) = smp ( N/2 * q + 1 : N/2 * q + N ) + d2 ;

while ( Pitch ( m1 ) ~= N/2 ) & ( m2 == 1 )
    m1=m1+1;
end

m2 = m2 + 1;
m1 = m1 + 1;

while ( m2 <= size(newfrms,1) ) & ( m1 < size ( Pitch , 2 ) )
    
    kk = newfrms( m2,:)*pinv(CBFBS);  %/ /(CBFBS);
    kk ( N : -1 : N-kc+1 ) = kk ( 2 : kc+1 );  % mirroring
    kk ( kc+2 : N-kc ) = 0;
    dd =  real ( ifft ( sqrt(kk) , N ) ) ;
    d = dd .* hamming (  N )' ;  %for deletting the effect of needle picks at first and end of frames 
%     d2=d;
    d2 = [d(N/2+1:N),d(1:N/2)]; 
    PP = Pitch ( m1 );
    
    while ( Pitch ( m1 ) ~= N/2 ) & ( m1 < size ( Pitch , 2 ) )
        smp ( N/2 * q + PP : N/2 * q + N + PP  -1) = smp ( N/2 * q + PP: N/2 * q + N + PP  -1) + d2 ;
        m1 = m1 + 1 ;
        PP = PP + Pitch ( m1 ) -1 ;
    end
    
    q = q + 1;    
    %     smp ( 512 * q  : 512 * q + 1024 - 1 ) = smp ( 512 * q  : 512* q + 1024 - 1 ) + d2 ;
    m2 = m2 + 1;
    m1 = m1 + 1;
    
end

wavwrite ( smp, fs , 'C:\Users\Admin\BME\project arshad\data\Wave\W4100_2.WAV' ) ;
% g2 = wavread ( 'D:/data/s' );
sound ( smp , fs );  
