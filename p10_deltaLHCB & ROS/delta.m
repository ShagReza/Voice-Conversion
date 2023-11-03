clc; close all; clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ROS signal dovom manand signal aval mishaad
mean_CB=load('C:\Users\Admin\BME\project arshad\programs\p5_mean & var extraction\mean_CB.mat');
mean_CB=mean_CB.mean_CB;
var_CB=load('C:\Users\Admin\BME\project arshad\programs\p5_mean & var extraction\var_CB.mat');
var_CB=var_CB.var_CB;

num_session=5;
num_subject=102;
marja=100;
z11=load(['C:\Users\Admin\BME\project arshad\data\LABLES\Zt',num2str(num_session),num2str(marja)]);
z1=z11.Z;%z1=label1
z22=load(['C:\Users\Admin\BME\project arshad\data\LABLES\Zt',num2str(num_session),num2str(num_subject)]);
z2=z22.Z;%z2=label2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shl=2; %entekhabe beyne signal1 o 2
while shl~=0
    N=512;%1024
    if shl==2
        [s1,fs]=wavread (['C:\Users\Admin\BME\project arshad\data\WAV\W',num2str(num_session),num2str(marja),'.WAV']);
    else
        filename = ['C:\Users\Admin\BME\project arshad\data\WAV\W',num2str(num_session),num2str(num_subject),'.wav'];
        [s1,fs] = wavread(filename);
    end
    % Pitch Extraction -------------------------------------------
    m = 1; Pitch=[]; PP=0;
    while m < length(s1) - N + 1,
        frm = s1( m : m + N - 1 );
        a = real ( ifft ( abs ( fft ( frm , N ) ) ) );
        [cc,bb] = max( a ( 50 : 256 ) ) ; %(100:512)
        pitch = bb + 49 ;%99
        m = m + pitch;
        Pitch=[Pitch pitch];
    end
    Pitch4=Pitch;
    %%%%%eslahe pitch vase sorat kond kardan%%%%%%%%%%%%%%%%%%%%%
    m = 1; PP=0; PITCH=[];
    while m < length(s1) - N + 1,
        frm = s1( m : m + N - 1 );
        a = real ( ifft ( abs ( fft ( frm , N ) ) ) );
        [cc,bb] = max ( a ( 50 : 256 ) ) ; %(100:512)
        pitch = bb + 49 ; %99
        m = m + 256; %m=m+512
        PITCH=[PITCH pitch];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m=1; i=1;
    while i <= length ( Pitch )
        m = m + Pitch ( i )  ;
        if m >= 256
            for j = length ( Pitch ) + 1 : -1 : i +1
                Pitch ( j ) = Pitch ( j - 1 ) ;
            end
            Pitch ( i ) = 256;
            Pitch ( i + 1 ) = m - 256 ;
            m = 0;
        end
        i = i + 1;
    end
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
    CBw1=CBw;
    cb11=cb(:,1:18).^2;
    cb1=sum(cb11);
    cbcepst11= cbcepst(:,1:18).^2;
    cbcepst1=sum(cbcepst11);

    if shl==2
        CB1=CBw;
        Pitch1=Pitch;
        CBFBS1l=CBFBS;
        clear CBw;
        clear Pitch;
        clear  frmfft;
    else
        CB2=CBw;
        Pitch2=Pitch;
        CBFBS2l=CBFBS;
    end
    shl=shl-1;
end %while shl~=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hazfe miyangin o variyans az LHCB:
mean_CB1=mean_CB(marja,:);
var_CB1=var_CB(marja,:);
CB1(:,1:18) = CB1(:,1:18) - repmat(mean_CB1 , size(CB1 , 1) ,1) ;
CB1(:,1:18) = CB1(:,1:18) ./ repmat(var_CB1 , size(CB1 , 1) ,1) ;
mean_CB2=mean_CB(num_subject,:);
var_CB2=var_CB(num_subject,:);
CB2(:,1:18) = CB2(:,1:18) - repmat(mean_CB2 , size(CB2 , 1) ,1) ;
CB2(:,1:18) = CB2(:,1:18) ./ repmat(var_CB2, size(CB2 , 1) ,1) ;

%extracting delte LHCB
size1=size(CB1);
size2=size(CB2);
delta2=[];
delta1=[];
for i=1:(size1(1)-1)
    delta1=[delta1,abs(mean(CB1(i,:)-CB1(i+1,:)))];
end
for i=1:(size2(1)-1)
    delta2=[delta2,abs(mean(CB2(i,:)-CB2(i+1,:)))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
index2=[]; index1=[];
DELTA1=delta1; DELTA2=delta2;
%mohasebeye maghadire maximomha va mahale anha dar delta baznamii sub1
i=1;
while i<45 %meghdare 45 ra ba tavajoh be tedad avahaye jomle moshakhas mikonim
    [a,b]=max(delta1);
    index1=[index1,b];
    delta1(b-2:b+2)=-1;
    i=i+1;
end
i=1;
while i<45
    [a,b]=max(delta2);
    index2=[index2,b];
    delta2(b-2:b+2)=-1;
    i=i+1;
end
index1=[index1,1,length(CB1)];
index2=[index2,1,length(CB2)];
ind1=sort(index1);
ind2=sort(index2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j1=1; j2=1;
npitch=[];
CB_Z2=[];
for ii=1:46
    len1=ind1(ii+1)-ind1(ii);
    len2=ind2(ii+1)-ind2(ii);
    %%%genetic algoritm for converting  to  :
    clear CB_label1 CB_label2
    CB_label1=CB1(j1:j1+len1-1,:);
    CB_label2=CB2(j2:j2+len2+1,:);
    j1=j1+1; j2=j2+1;
    %tolid jamiyate avaliye:
    clear CB_new_total k CB_new label_new CB_new_total label_new_half label_motivation label_new_total LABEL LABEL_mot LABEL_half

    if len1>len2
        for i=1:100
            label_new=1:len2;
            L=randint(1,len1-len2,[1,len2]);
            label_new=[label_new,L];
            LABEL(i,1:len1)=label_new;
            label_new=sort(label_new); %baraye soodi bodane dad
            label_new_total(i,1:len1)=label_new;
        end
        %%%
        for h=1:10
            clear CB_new_total ;
            %CB_new_total:
            for i=1:100
                for j=1:len1
                    k=label_new_total(i,j);
                    CB_new(j,1:175)=CB_label2(k,1:175);
                end
                CB_new_total((i-1)*len1+1:i*len1,1:175)=CB_new;
            end
            %taabe barazandegi:
            clear error
            for i=1:100
                CB_new= CB_new_total((i-1)*len1+1:i*len1,1:175);
                error(i)=sum(sum((CB_new-CB_label1).*(CB_new-CB_label1)));
            end
            clear E
            E=sort(error); j=0;
            %entekhabe barazandeha:
            clear label_new_half,clear LABEL_half
            for i=1:100
                if error(i)<=E(50)
                    j=j+1;
                    label_new_half(j,:)=label_new_total(i,:);
                    LABEL_half(j,:)=LABEL(i,:);
                end
            end
            %motivation:
            clear label_motivation,clear LABEL_mot
            label_motivation(1:50,1:len1)=LABEL_half(1:50,1:len1);
            X=min(len1,len2);
            for i=1:50
                k1=randint(1,1,[len2+1,len1]);
                k2=randint(1,1,[1,len2]);
                label_motivation(i,k1)=k2;
                LABEL_mot(i,:)=label_motivation(i,:);
                label_motivation(i,:)=sort(label_motivation(i,:));
            end
            clear label_new_total,clear LABEL
            label_new_total=[label_new_half;label_motivation];
            LABEL=[LABEL_half;LABEL_mot];
        end %for h=1:20
        [a,b]=min(error);
        label=label_new_total(b,:)
    end %if len1>len2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if len1<len2
        for i=1:100
            label_new=randint(1,len1,[1,len2]);
            label_new=sort(label_new); %baraye soodi bodane dad
            label_new_total(i,1:len1)=label_new;
        end
        for h=1:40
            clear CB_new_total ;
            %CB_new_total:
            for i=1:100
                for j=1:len1
                    k=label_new_total(i,j);
                    CB_new(j,1:175)=CB_label2(k,1:175);
                end
                CB_new_total((i-1)*len1+1:i*len1,1:175)=CB_new;
            end
            %taabe barazandegi:
            clear error
            for i=1:100
                CB_new= CB_new_total((i-1)*len1+1:i*len1,1:175);
                error(i)=sum(sum((CB_new-CB_label1).*(CB_new-CB_label1)));
            end
            clear E
            E=sort(error); j=0;
            %entekhabe barazandeha:
            clear label_new_half
            for i=1:100
                if error(i)<=E(50)
                    j=j+1;
                    label_new_half(j,:)=label_new_total(i,:);
                end
            end
            %tolid nasle bad:
            %crossover:
            m=floor(len1/2);
            clear label_crossover
            label_crossover(1:25,1:len1)=0;
            for i=1:25
                label_crossover(i,:)=[label_new_half(i,1:m) label_new_half(i+1,m+1:len1)];
            end
            %motivation:
            clear label_motivation
            label_motivation(1:25,1:len1)=label_new_total(1:25,1:len1);
            X=min(len1,len2);
            for i=1:25
                k1=randint(1,1,[1,len1]);
                k2=randint(1,1,[1,len2]);
                label_motivation(i,k1)=k2;
                label_motivation(i,:)=sort(label_motivation(i,:));

                k1=randint(1,1,[1,len1]);
                k2=randint(1,1,[1,len2]);
                label_motivation(i,k1)=k2;
                label_motivation(i,:)=sort(label_motivation(i,:));

                k1=randint(1,1,[1,len1]);
                k2=randint(1,1,[1,len2]);
                label_motivation(i,k1)=k2;
                label_motivation(i,:)=sort(label_motivation(i,:));
            end
            clear label_new_total
            label_new_total=[label_new_half; label_crossover; label_motivation];
        end %for h=1:20
        [a,b]=min(error);
        label=sort(label_new_total(b,:))
    end %if len1<=len2
    if len1==len2
        label=1:len1
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PITCH_len2=PITCH(ind2(ii):ind2(ii+1));
    PITCH_label=PITCH_len2(1);
    for i=2:len1
        j=label(i);
        if label(i)==label(i-1) && j~=len2
            PITCH_label(i)=floor([PITCH_len2(j+1)+PITCH_len2(j)]/2);
        else
            PITCH_label(i)=PITCH_len2(j);
        end
    end
    npitch=[npitch PITCH_label];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CB_new=[];
    for x=1:len1
        k=label(x);
        CB_new(x,1:175)=CB_label2(k,1:175);
    end
    CB_Z2=[CB_Z2;CB_new];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mm=0;pitch_new_Z2=[];jj=0;
for i=1:length(npitch)
    while mm<256 %512
        mm=mm+npitch(i);
        if mm<256 %512
            jj=jj+1;
            pitch_new_Z2(jj)=npitch(i);
        else
            pitch_new_Z2(jj+1)=256; %512
            pitch_new_Z2(jj+2)=mm-256; %512
        end
    end
    mm=mm-256;%512
    jj=jj+2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CB_Z2 to wave:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ezafe kardane miyangin o variyans az LHCB:
CB_Z2(:,1:18) = CB_Z2(:,1:18) .* repmat(var_CB2, size(CB1 , 1)-1 ,1) ;
CB_Z2(:,1:18) = CB_Z2(:,1:18) + repmat(mean_CB2 , size(CB1 , 1)-1,1) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pitch=pitch_new_Z2; %!!!
m1=1; N2=175; m2=1; CBw=CB_Z2;
newfrms = exp ( CBw ) - 0.01 ;
z=size(newfrms,1)*512-(size(newfrms,1)-1)*256;
smp=zeros(1,z);
PP=0;  q=0;

kk = newfrms( m2,:) / CBFBS; %!!!
kk ( 512 : -1 : 339 ) = kk ( 2 : 175 );  % mirrorong
kk ( 176 : 340 ) = 0;
dd =  real ( ifft ( sqrt(kk) , 512 ) ) ;
d = dd .* hamming(512)' ;  %for deletting the effect of needle picks at first and end of frames
d2 = [d(257:512),d(1:256)];
smp ( 256 * q + 1 : 256 * q + 512 ) = smp ( 256 * q + 1 : 256 * q + 512 ) + d2 ;

while ( Pitch ( m1 ) ~= 256 ) & ( m2 == 1 )
    m1=m1+1;
end

m2 = m2 + 1; m1 = m1 + 1;

while ( m2 <= size(newfrms,1) ) & ( m1 < size ( Pitch , 2 ) )

    %kk = newfrms( m2,:) / CBFBS; %!!!
    kk = newfrms( m2,:)*pinv(CBFBS);  %/(CBFBS);
    kk ( 512 : -1 : 339 ) = kk ( 2 : 175 );  % mirrorong
    kk ( 176 : 340 ) = 0;
    dd =  real ( ifft ( sqrt(kk) , 512 ) ) ;
    d = dd .* hamming(512)' ;  %for deletting the effect of needle picks at first and end of frames
    d2 = [d(257:512),d(1:256)];

    PP = Pitch ( m1 );

    while ( Pitch ( m1 ) ~= 256 ) & ( m1 < size ( Pitch , 2 ) )
        smp ( 256 * q + PP : 256 * q + 512 + PP  -1) = smp ( 256 * q + PP : 256 * q + 512 + PP  -1) + d2 ;
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
