%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% In the name of God %%%%%%%
%%%%%%% shaghayegh reza %%%%%%%%%
%%%%%%%% project arshad:ROS changing %%%%%%%%%
%%%%%%%%%% 27 farvardin 87 %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; close all; clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ROS signal dovom manand signal aval mishaad
mean_CB=load('Z:\project arshad\data\mean_CB.mat');
mean_CB=mean_CB.mean_CB;
var_CB=load('Z:\project arshad\data\var_CB.mat');
var_CB=var_CB.var_CB;
for num_session=4:5
    for num_subject=100:200
        marja=113;
        z11=load(['Z:\project arshad\data\LABLES\Zt',num2str(num_session),num2str(marja)]);
        z1=z11.Z;%z1=label1
        z22=load(['Z:\project arshad\data\LABLES\Zt',num2str(num_session),num2str(num_subject)]);
        z2=z22.Z;%z2=label2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        shl=2; %entekhabe beyne signal1 o 2
        while shl~=0
            N=1024; n=5; K=105;%184
            if shl==2
                [s1,fs]=wavread ([ 'Z:\project arshad\data\WAV\W',num2str(num_session),num2str(marja),'.WAV']);
            else
                filename = ['Z:\project arshad\data\Wav\W',num2str(num_session),num2str(num_subject),'.wav'];
                [s1,fs] = wavread(filename);
            end
            % Pitch Extraction -------------------------------------------
            m = 1; Pitch=[]; PP=0;
            while m < length(s1) - N + 1,
                frm = s1( m : m + N - 1 );
                a = real ( ifft ( abs ( fft ( frm , N ) ) ) );
                [cc,bb] = max ( a ( 100 : 512 ) ) ;
                pitch = bb + 99 ; m = m + pitch;
                Pitch=[Pitch pitch];
            end
            Pitch4=Pitch;
            %%%%%eslahe pitch vase sorat kond kardan%%%%%%%%%%%%%%%%%%%%%
            m = 1; PP=0; PITCH=[];
            while m < length(s1) - N + 1,
                frm = s1( m : m + N - 1 );
                a = real ( ifft ( abs ( fft ( frm , N ) ) ) );
                [cc,bb] = max ( a ( 100 : 512 ) ) ;
                pitch = bb + 99 ; m = m + 512;
                PITCH=[PITCH pitch];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            m=1; i=1;
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
            %-feature extraction--------------------------------------------------
            N = 1024 ;		            % Frame length
            stp = 512;		            % Frame rate
            L = 18 ;		            % Number of LHCB parameters
            Q = 18 ;		            % Number of Filters
            fc = 7500 ;	                % Final frequency is 7500 Hz
            fs = 44100 ;            	% Sampling frequency
            kc = round( fc * N / fs );	% Equivallent of 7500 Hz in DFT domain
            f = fs / N * (0:kc);		% Frequencies samples in hertz
            win = hamming(N);	        % Temporal window
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
                frmfftsqr = real ( frmfft .* conj ( frmfft ) )';
                %energy of filter
                cb = frmfftsqr * CBFBS ;
                cbfinal1=[cbfinal1 cb];
                cbcepst = log ( 0.01 + cb ) ;
                cbfinal2=[cbfinal2 cbcepst];
                %after log
                CBw = [CBw ; cbcepst] ;
                m = m + 512 ;  n = n + 1;
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
%         mean_CB1=mean_CB(marja,:);
%         var_CB1=var_CB(marja,:);
%         CB1(:,1:18) = CB1(:,1:18) - repmat(mean_CB1 , size(CB1 , 1) ,1) ;
%         CB1(:,1:18) = CB1(:,1:18) ./ repmat(var_CB1 , size(CB1 , 1) ,1) ;
%         mean_CB2=mean_CB(num_subject,:);
%         var_CB2=var_CB(num_subject,:);
%         CB2(:,1:18) = CB2(:,1:18) - repmat(mean_CB2 , size(CB2 , 1) ,1) ;
%         CB2(:,1:18) = CB2(:,1:18) ./ repmat(var_CB2, size(CB2 , 1) ,1) ;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %amaliate hamahang sazie 2 seda%%%%%%%%%%%%%
        Pitch1=[Pitch1 N/2]; Pitch2=[Pitch2 N/2]; Pitch3=Pitch2;
        op=[]; ol=[]; fd=0; cu=0; test=0; dsl=length(z1);
        ad=1; bsl=1;par2=zeros(dsl,175); par1=zeros(dsl,175);
        par1=CB1; npitch1=Pitch1; co=0; cue=0; AAA=0;

        % tashkhise entehaye har frame dar matris Pitch
        while co~=length(z2)
            cue=cue+1;
            if Pitch2(1,cue)==N/2
                co=co+1;
                cp(1,co)=cue;
            end
        end
        %%%%%%%%%%
        cpp=cp+1;
        index_pitch=[1 cpp];
        Pitch3=[];
        %%%%%%%%%%
        mm=0;npitch=[];jj=0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % dar in halghe toole avahaye motanazer tashkhis dade shode
        len_z1=length(z1); %lenght of label of signal1
        len_z2=length(z2);
        z1_counter=1; z2_counter=1; %z1_counter is the first & z1_count t40=0;e end of label
        z1_count=1; z2_count=1;  CB_Z2=[]; W=0; len1=0; len2=0; END=1; Alarm40=0; Q=0;

        while z1_count<(len_z1) && END==1
            %% counting label1 of z1________________________________________
            z1_counter=z1_count; label1_z1=z1(z1_counter,1);   count_z1_label1=0; %lenght of label1 of z1
            while z1(z1_count,1)==label1_z1 && z1_count~=(len_z1)
                count_z1_label1=count_z1_label1+1;
                z1_count=z1_count+1;
            end
            %% counting label1 of z2
            z2_counter=z2_count; label1_z2=z2(z2_counter,1);    count_z2_label1=0; %lenght of label1 of z1
            while z2(z2_count,1)==label1_z2 && z2_count~=(len_z2)
                count_z2_label1=count_z2_label1+1;
                z2_count=z2_count+1;
            end
            %%_____________________________________________________________
            if z1_count==len_z1 || z2_count==len_z2
                len2=length(z2_counter:len_z2);
                len1=length(z1_counter:len_z1);
                END=0;
            else
                len1=count_z1_label1;
                len2=count_z2_label1;
            end

            if label1_z1~=label1_z2 && END==1
                %%_______________reading the next label:___________________
                %% counting label2 of z1
                z1_count2=z1_count; label2_z1=z1(z1_count2,1);  count_z1_label2=0;
                while z1(z1_count2,1)==label2_z1 && z1_count2~=(len_z1)
                    count_z1_label2=count_z1_label2+1;
                    z1_count2=z1_count2+1;
                end
                %% counting label2 of z2
                z2_count2=z2_count;  label2_z2=z2(z2_count,1);  count_z2_label2=0;
                while z2(z2_count2,1)==label2_z2 &&  z2_count2~=(len_z2)
                    count_z2_label2=count_z2_label2+1;
                    z2_count2=z2_count2+1;
                end
                %%_________________________________________________________
                if z1_count2==len_z1 || z2_count2==len_z2
                    len2=length(z2_counter:len_z2);
                    len1=length(z1_counter:len_z1);
                    END=0;
                else
                    len1=count_z1_label1; len2=count_z2_label1;
                end
                %%%%%%%%%%
                if label2_z1~=label2_z2 && END==1
                    if label2_z1==label1_z2 && label1_z1==label2_z2
                        len1=count_z1_label2+count_z1_label1;
                        len2=count_z2_label2+count_z2_label1;
                        z1_count=z1_count2; z2_count=z2_count2;
                    elseif label2_z1==label1_z2
                        if label1_z1==40 %!!!!!!!!!!!
                            len1=count_z1_label1; len2=len1; z2_count=z2_counter; Alarm40=1;
                        else
                            len1=count_z1_label2+count_z1_label1; len2=count_z2_label1;
                            z1_count=z1_count2;
                        end
                    elseif label1_z1==label2_z2
                        if label1_z2 == 40
                            len1=count_z1_label1; len2=count_z2_label2;
                            W=W+count_z2_label1; z2_counter=z2_count; z2_count=z2_count2;
                        else
                            len1=count_z1_label1; len2=count_z2_label2+count_z2_label1;
                            z2_counter=z2_count; z2_count=z2_count2;
                        end
                    else
                        %___________reading the third label________________
                        %% counting label3 of z1
                        z1_count3=z1_count2; label3_z1=z1(z1_count2,1);count_z1_label3=0;
                        while z1(z1_count3,1)==label3_z1 && z1_count3<(len_z1)
                            count_z1_label3=count_z1_label3+1;
                            z1_count3=z1_count3+1;
                        end
                        %% counting label2 of z2
                        z2_count3=z2_count2; label3_z2=z2(z2_count2,1);  count_z2_label3=0;
                        while z2(z2_count3,1)==label3_z2 && z2_count3<(len_z2)
                            count_z2_label3=count_z2_label3+1;
                            z2_count3=z2_count3+1;
                        end
                        %___________________________________________________________
                        if label3_z1==label3_z2
                            if label1_z1==40
                                len1=count_z1_label1; len2=len1; z2_count=z2_counter; Alarm40=1;
                            elseif label1_z2==40 %faghat 40 hazf va dobare khande shavad
                                Alarm40=1;
                                W=W+count_z2_label1;
                                z1_count=z1_counter;
                            else
                                len1=count_z1_label1+count_z1_label2; len2=count_z2_label1+count_z2_label2;
                                z1_count=z1_count2; z2_count=z2_count2;
                            end
                        elseif label3_z1==label1_z2%ok
                            len1=count_z1_label1+count_z1_label2; len2=len1;
                            Q=1; z1_count=z1_count2; z2_count=z2_counter;
                        elseif label3_z1==label2_z2
                            len1=count_z1_label1+count_z1_label2; len2=count_z2_label1; z1_count=z1_count2;
                        elseif label1_z1==label3_z2%ok
                            len1=count_z1_label1; len2=count_z2_label3; z2_counter=z2_count2; z2_count=z2_count3;
                            W=W+count_z2_label1+count_z2_label2;
                        elseif label2_z1==label3_z2
                            len1=count_z1_label1; len2=count_z2_label1+count_z2_label2; z2_count=z2_count2;
                        else %tasavi vojood nadasht fagat 2 label aval tabdil shavand
                            len1=count_z1_label1; len2=count_z2_label1;
                        end
                    end
                end % if label2_z1==label2_z2
            end %if label1_z1~=label1_z2
            len1,len2
            %%%genetic algoritm for converting  to  :
            clear CB_label1 CB_label2
            CB_label1=CB1(z1_counter:z1_counter+len1-1,:);
            CB_label2=CB2(z2_counter:z2_counter+len2-1,:);
            %tolid jamiyate avaliye:
            clear CB_new_total k CB_new label_new CB_new_total label_new_half label_motivation label_new_total LABEL LABEL_mot LABEL_half

            if len1>len2 && Alarm40==0 && Q==0
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
            if len1<len2 && Alarm40==0 && Q==0
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
                label=label_new_total(b,:)
            end %if len1<=len2
            if len1==len2
                label=1:len1
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if   Alarm40==0 && Q==0
                W=W+len2;
            elseif (Alarm40==1 && label1_z1==40) || (Q==1)
                ooo='oooALARM40!ooo'
                CB_label2=CB1(z1_counter:z1_counter+len1-1,:);
                XX=PITCH(1:len1);
            else
                XX=[];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if  Alarm40==0 && Q==0
                PITCH_len2=PITCH(z2_counter:z2_counter+len2-1);
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
            else
                npitch=[npitch XX];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if   (Alarm40==0 && Q==0) || (Alarm40==1 && label1_z1==40) || (Q==1)
                CB_new=[];
                for x=1:len1
                    k=label(x);
                    CB_new(x,1:175)=CB_label2(k,1:175);
                end
                CB_Z2=[CB_Z2;CB_new];
            end
            Alarm40=0; Q=0;
        end %while z1_counter<(len_z1-1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mm=0;pitch_new_Z2=[];jj=0;
        for i=1:length(npitch)
            while mm<512
                mm=mm+npitch(i);
                if mm<512
                    jj=jj+1;
                    pitch_new_Z2(jj)=npitch(i);
                else
                    pitch_new_Z2(jj+1)=512;
                    pitch_new_Z2(jj+2)=mm-512;
                end
            end
            mm=mm-512;
            jj=jj+2;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% CB_Z2 to wave:
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %ezafe kardane miyangin o variyans az LHCB:
        %CB_Z2(:,1:18) = CB_Z2(:,1:18) .* repmat(var_CB2, size(CB1 , 1) ,1) ;
        %CB_Z2(:,1:18) = CB_Z2(:,1:18) + repmat(mean_CB2 , size(CB1 , 1) ,1) ;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Pitch=pitch_new_Z2; %!!!
        m1=1; N2=175; m2=1; CBw=CB_Z2;
        newfrms = exp ( CBw ) - 0.01 ;
        z=size(newfrms,1)*1024-(size(newfrms,1)-1)*512;
        smp=zeros(1,z);
        PP=0;  q=0;

         %kk = newfrms( m2,:) / CBFBS; %!!!
            kk = newfrms( m2,:)*pinv(CBFBS);  %/(CBFBS);
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

            %kk = newfrms( m2,:) / CBFBS; %!!!
            kk = newfrms( m2,:)*pinv(CBFBS);  %/(CBFBS);
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
        %%%%%%%%mohasebe energy:%%%%%%%%%%
        energy_marja=sum(s1.*s1);
        energy_smp=sum(smp.*smp);
        smp=smp*[energy_marja/energy_smp];
        %smp=smp*5;
        wavwrite ( smp, fs ,['Z:\project arshad\data\equal ROS waves\W',num2str(num_session),num2str(num_subject),'_to_',num2str(marja),'.wav']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %makhlot kardane do soot:
        [s1,fs]=wavread ([ 'Z:\project arshad\data\WAV\W',num2str(num_session),num2str(marja),'.WAV']);
        S=0;
        s1=s1(1:length(smp));
        S=s1+smp';
        wavwrite ( S', fs ,['Z:\project arshad\data\equal ROS waves\W',num2str(num_session),num2str(num_subject),'_makhloot_',num2str(marja),'.wav']);
    end % for num_subject=100:200
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

