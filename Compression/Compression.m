%% LOSSY COMPRESSION-DECOMPRESSION USNIG DISCRETE COSINE TRANSFORM CALCULATOR BY FFT.
clear;
clc;
[filename,pathname] = uigetfile('*.*','Select grayscale Image');
filewithpath = strcat(pathname,filename);
I = imread(filename);
n = 8;
            m = 7;
            
            
            N=8;                        % Block size for which DCT is C/omputed.
            M=8;
                % ??c h?nh ??u v?o.
            
            
            if length(size(I))>2
                I=rgb2gray(I);
            end
            ori= I;
            I_dim_old=size(I);          % ?? d?i c?a b?c ?nh.
            I_dim_new = I_dim_old ;
            
            
            
            % Ch?nh s?a ?? d?i ?nh t?o bi?n ?o
            orinal_x = log2(I_dim_old(1));
            orinal_y = log2(I_dim_old(2));
            while mod(log2(I_dim_new(1)),1) ~= 0
                I_dim_new(1) = 2^(orinal_x - mod(orinal_x,1) + 1) ;
                for x = I_dim_old(1):I_dim_new(1)
                    I(x,:) = I(I_dim_old(1),:);
                end
            end
            while mod(log2(I_dim_new(2)),1) ~= 0 
                I_dim_new(2) = 2^(orinal_y - mod( orinal_y,1) + 1);
                for y = I_dim_old(2):I_dim_new(2)
                    I(:,y) = I(:,I_dim_old(2));
                end
            end
            while I_dim_new(1) < I_dim_new(2)
                c = I_dim_new(1);
                I_dim_new(1) = I_dim_new(2);
                for x = c :I_dim_new(2)
                    I(x,:) = I(I_dim_old(1),:);
                end
            end
            while I_dim_new(1) > I_dim_new(2) 
                c = I_dim_new(2);
                for y = I_dim_new(2):I_dim_new(1)
                    I(:,y) = I(:,I_dim_old(2));
                end
            end
            
            % L?u gi? tr? h?nh m?i v?o bi?n I_dim.
            I_dim = I_dim_new ;
            mask = zeros(8);
            I_Trsfrm.block=zeros(N,M);  % Kh?i t?o ma tr?n Block 8x8.
            
            
            
            k = 20;
            
            
            
            
            
            
            Norm_Mat=[16 11 10 16 24 40 51 61       % Normalization matrix (8 X 8) used to Normalize the DCT Matrix.
                      12 12 14 19 26 58 60 55
                      14 13 16 24 40 57 69 56
                      14 17 22 29 51 87 80 62
                      18 22 37 56 68 109 103 77
                      24 35 55 64 81 104 113 92
                      49 64 78 87 103 121 120 101
                      72 92 95 98 112 100 103 99];
             
                  
             Norm_Mat=Norm_Mat*k;     
            save('Initial.txt','I');
            
            
        
            imwrite(I,'orani.jpg');
            
            
            
            
            
            %% PART-1: COMPRESSION 
            
            
            for a=1:I_dim(1)/N
                for b=1:I_dim(2)/M
                    for i=1:N
                        for j=1:M
                            mask(i,j) = I(N*(a-1)+i,M*(b-1)+j);
                        end
                    end
                    I_Trsfrm(a,b).block =  dct2d(mask);
                    % Quantization Table
                    I_Trsfrm(a,b).block = round(I_Trsfrm(a,b).block./Norm_Mat);
                end
            end
            
            %Thu?t to?n zig zag cho ma tr?n 8x8 .
            for a=1:I_dim(1)/N
                for b=1:I_dim(2)/M
                    I_zigzag(a,b).block=zeros(1,0);
                    freq_sum=2:(N+M);
                    counter=1;
                    for i=1:length(freq_sum)
                        if i<=((length(freq_sum)+1)/2)
                            if rem(i,2)~=0
                                x_indices=counter:freq_sum(i)-counter;
                            else
                                x_indices=freq_sum(i)-counter:-1:counter;
                            end
                                index_len=length(x_indices);
                                y_indices=x_indices(index_len:-1:1); % Creating reverse of the array as "y_indices".
                                for p=1:index_len
                                    if I_Trsfrm(a,b).block(x_indices(p),y_indices(p))<0
                                        bin_eq=dec2bin(bitxor(2^n-1,abs(I_Trsfrm(a,b).block(x_indices(p),y_indices(p)))),n);
                                    else
                                        bin_eq=dec2bin(I_Trsfrm(a,b).block(x_indices(p),y_indices(p)),n);
                                    end
                                    I_zigzag(a,b).block=[I_zigzag(a,b).block,bin_eq(1:m)];
                                end
                        else
                            counter=counter+1;
                            if rem(i,2)~=0
                                x_indices=counter:freq_sum(i)-counter;
                            else
                                x_indices=freq_sum(i)-counter:-1:counter;
                            end
                                index_len=length(x_indices);
                                y_indices=x_indices(index_len:-1:1); % Creating reverse of the array as "y_indices".
                                for p=1:index_len
                                    if I_Trsfrm(a,b).block(x_indices(p),y_indices(p))<0
                                        bin_eq=dec2bin(bitxor(2^n-1,abs(I_Trsfrm(a,b).block(x_indices(p),y_indices(p)))),n);
                                    else
                                        bin_eq=dec2bin(I_Trsfrm(a,b).block(x_indices(p),y_indices(p)),n);
                                    end
                                    I_zigzag(a,b).block=[I_zigzag(a,b).block,bin_eq(1:m)];
                                end
                        end
                    end
                end
            end
            
            % Clearing unused variables from Memory space
            clear I_Trsfrm prod; 
            clear x_indices y_indices counter;
            
            % Run-Length Encoding.
            for a=1:I_dim(1)/N
                for b=1:I_dim(2)/M
                    
                    % T?nh to?n Count cho c?c k? hi?u t??ng ?ng v? 
                    % l?u ch?ng v?o "I_run" structure.
                    count=0;
                    run=zeros(1,0);
                    sym=I_zigzag(a,b).block(1);
                    j=1;
                    block_len=length(I_zigzag(a,b).block);
                    for i=1:block_len
                        if I_zigzag(a,b).block(i)== sym
                            count=count+1;
                        else
                            run.count(j)=count;
                            run.sym(j)=sym;
                            j=j+1;
                            sym=I_zigzag(a,b).block(i);
                            count=1;
                        end
                        if i==block_len
                            run.count(j)=count;
                            run.sym(j)=sym;
                        end
                    end 
                    
                    %T?nh to?n c??ng ?? runlength code cho c?c gi? tr? ??m 
                    dim=length(run.count);  % S? l??ng k? hi?u ???c m? h?a .
                    maxvalue=max(run.count);  % T?m gi? tr? ??m l?n nh?t trong m?ng.
                    codelength=log2(maxvalue)+1;
                    codelength=floor(codelength);
                    
                    % M? h?a c?c gi? tr? ??m c?ng v?i c?c k? hi?u .
                    I_runcode(a,b).code=zeros(1,0);
                    for i=1:dim
                        I_runcode(a,b).code=[I_runcode(a,b).code,dec2bin(run.count(i),codelength),run.sym(i)];
                    end
                end
            end
            % L?u m? n?n .
            Rate = zeros(1,0);
            for a= 1:I_dim(1)/N
                for b = 1:I_dim(2)/M
                    Rate = [Rate,I_runcode(a,b).code]; 
                end
            end
            Rate_test = char(Rate);
            save ('Compressed.txt',"Rate_test");
            
            
            
            
            fileID = fopen('myfile.txt','w'); 
            
            nbytes = fprintf(fileID,Rate);
            fclose(fileID);
            rate= 0;
            for a=1:I_dim(1)/N
                for b=1:I_dim(2)/M
                  rate = rate +  length(I_runcode(a,b).code)/8000;  
                
                end
            end
            
%% PART-2: DECOMPRESSION TECHNIQUE.

% Run-length Decoding .
for a=1:I_dim(1)/N
    for b=1:I_dim(2)/M
        enc_str=I_runcode(a,b).code;
        
        % T?nh ?? d?i c?a chu?i ???c m? h?a.
        enc_len=length(enc_str);
        
        % Since Max. Count is unknown at the receiver, Number of bits used for each 
        % count value is unknown and hence cannot be decoded directly. Number of bits 
        % used for each count can be found out by trial and error method for all 
        % the possible lengths => factors of encoded string length.

        % T?nh to?n the non-trivial factors of the "enc_len" (length of encoded string) i.e., factors other than 1 & itself.
        factors_mat=zeros(1,0);
        if enc_len<=(n+1)
            realfact=enc_len;
        else
            for i=2:enc_len-2       % "enc_len-1" Lu?n lu?n kh?ng ph?i l? ??c s? "enc_len".
                if(rem(enc_len,i)==0)
                    factors_mat=[factors_mat,i];
                end
            end

            % Th? v? l?i ?? t?m ra gi? tr? ??m ch?nh x?c.
            for i=1:length(factors_mat)
                flagcntr=0;
                temp_dim=enc_len/factors_mat(i);
                for j=1:temp_dim
                    if strcmp(enc_str(1+(j-1)*factors_mat(i):j*factors_mat(i)),dec2bin(0,factors_mat(i)))==0
                        if j==1
                            flagcntr=flagcntr+1;
                        else
                            if enc_str((j-1)*factors_mat(i))~=enc_str(j*factors_mat(i))
                                flagcntr=flagcntr+1;
                            else
                                break;
                            end
                        end
                    else
                        break;
                    end
                end
                if flagcntr==temp_dim
                    realfact=factors_mat(i);
                    break;
                end
            end
        end
        
        % Clearing unused variables from Memory space
        % clear factors_mat flagcntr j 

% T?m chu?i m? h?a ra c?c gi? tr? ??m c?a c?c k? hi?u t??ng ?ng
        dec_str=zeros(1,0);
        temp_dim=enc_len/realfact;
        for i=1:temp_dim
            count_str=enc_str(1+(i-1)*realfact:(i*realfact)-1);
            countval=bin2dec(count_str);
            for j=1:countval
                dec_str=[dec_str,enc_str(i*realfact)];
            end
        end
        I_runcode(a,b).code=dec_str;
    end
end

% D?n d?p kh?ng gian b? nh?
% clear enc_str dec_str temp_dim realfact enc_len
% clear countval count_str

% T?i t?o c?c kh?i block theo ki?u zig-zag.
I_rec_Trnsfm.block=zeros(N,M);
for a=1:I_dim(1)/N
    for b=1:I_dim(2)/M
        bpp=length(I_zigzag(a,b).block)/(N*M);  % "bpp" is the bits-per-pixel in reconstruction of image.
        bpp_diff= n-bpp;
        freq_sum=2:(N+M);
        counter=1;
        c_indx=1;
        for i=1:length(freq_sum)
            if i<=((length(freq_sum)+1)/2)
                if rem(i,2)~=0
                    x_indices=counter:freq_sum(i)-counter;
                else
                    x_indices=freq_sum(i)-counter:-1:counter;
                end
                    index_len=length(x_indices);
                    y_indices=x_indices(index_len:-1:1); % Creating reverse of the array as "y_indices".
                    for p=1:index_len
                        decm_eq=bin2dec([I_runcode(a,b).code(1+m*(c_indx-1):m*c_indx),dec2bin(0,bpp_diff)]);
                        if decm_eq>(2^(n-1))-1
                            decm_eq=decm_eq-(2^n-1);
                        end
                        I_rec_Trnsfm(a,b).block(x_indices(p),y_indices(p))=decm_eq;
                       c_indx=c_indx+1;
                    end
            else
                counter=counter+1;
                if rem(i,2)~=0
                    x_indices=counter:freq_sum(i)-counter;
                else
                    x_indices=freq_sum(i)-counter:-1:counter;
                end
                    index_len=length(x_indices);
                    y_indices=x_indices(index_len:-1:1); % Creating reverse of the array as "y_indices".
                    for p=1:index_len
                        decm_eq=bin2dec([I_runcode(a,b).code(1+m*(c_indx-1):m*c_indx),dec2bin(0,bpp_diff)]);
                        if decm_eq>(2^(n-1))-1
                            decm_eq=decm_eq-(2^n-1);
                        end
                        I_rec_Trnsfm(a,b).block(x_indices(p),y_indices(p))=decm_eq;
                        c_indx=c_indx+1;
                    end
            end
        end
    end
end

% D?n d?p kh?ng gian b? nh?
 clear I_runcode x_indices y_indices
 clear c_indx freq_sum


 
% Chu?n h?a ma tr?n t?i t?o b?ng Lumimnance Table.
for a=1:I_dim(1)/N
    for b=1:I_dim(2)/M
        I_rec_Trnsfm(a,b).block=(I_rec_Trnsfm(a,b).block).*Norm_Mat;
    end
end


% IDCT by IFFT.


for a=1:I_dim(1)/N
    for b=1:I_dim(2)/M
        mask = idct2d(I_rec_Trnsfm(a,b).block);
        for i=1:N
            for j=1:M
               I_rec((a-1)*N+i,(b-1)*M+j) = mask(i,j); 
            end
        end
    end
end

% D?n d?p kh?ng gian b? nh?
 clear I_rec_Trnsfm


% In ra m?n h?nh b?c ?nh.

I_rec=I_rec/max(max(I_rec));
I_rec=im2uint8(I_rec);
% C?t v?ng bi?n ?o
I = imcrop(I,[0 0 I_dim_old(2) I_dim_old(1)]);
I_rec= imcrop(I_rec,[0 0 I_dim_old(2) I_dim_old(1)]);
figure,imshow(I_rec,[0,2^n-1]);
figure, imhist(I_rec);
imwrite(I_rec,'Compress.jpg');
psnr(I,I_rec);