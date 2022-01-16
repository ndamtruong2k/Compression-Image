function varargout = nenanh(varargin)
% NENANH MATLAB code for nenanh.fig
%      NENANH, by itself, creates a new NENANH or raises the existing
%      singleton*.
%
%      H = NENANH returns the handle to a new NENANH or the handle to
%      the existing singleton*.
%
%      NENANH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NENANH.M with the given input arguments.
%
%      NENANH('Property','Value',...) creates a new NENANH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before nenanh_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to nenanh_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help nenanh

% Last Modified by GUIDE v2.5 31-Mar-2021 16:16:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nenanh_OpeningFcn, ...
                   'gui_OutputFcn',  @nenanh_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before nenanh is made visible.
function nenanh_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nenanh (see VARARGIN)

% Choose default command line output for nenanh
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes nenanh wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = nenanh_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in OPEN BUTTON.



function pushbutton3_Callback(hObject, eventdata, handles)

[filename, pathname] = uigetfile('*.*', 'Select grayscale Image');
    if isequal(filename,0) || isequal(pathname,0)
       disp('User pressed cancel')
    else
       filename=strcat(pathname,filename);
       global I;
       I=imread(filename);
       axes(handles.axes1);
       imshow(I);
       handles.o=I;
      
       guidata(hObject, handles);
    end
    
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)






% --- Executes on button press in COMPRESSION BUTTON.

function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)


%% LOSSY COMPRESSION-DECOMPRESSION USNIG DISCRETE COSINE TRANSFORM CALCULATOR BY FFT.


clear;
clc;
global I;
n = 8;
m = 7;




N=8;                        % Block size for which DCT is C/omputed.
M=8;
    % ??c h?nh ??u v?o.


if length(size(I))>2
    I=rgb2gray(I);
end;
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
    I(I_dim_new(1):1:I_dim_new(2),:) = I(I_dim_old(1),:);
    for x = I_dim_new(1):I_dim_new(2)
        I(x,:) = I(I_dim_old(1),:);
    end
    I_dim_new(1) = I_dim_new(2);
end
while I_dim_new(1) > I_dim_new(2) 
    I(:,I_dim_new(2):1:I_dim_new(1)) = I(:,I_dim_old(2));
    for y = I_dim_new(2):I_dim_new(1)
        I(:,y) = I(:,I_dim_old(2));
    end
    I_dim_new(2) = I_dim_new(1) ;
end

% L?u gi? tr? h?nh m?i v?o bi?n I_dim.
I_dim = I_dim_new ;
mask = zeros(8);
I_Trsfrm.block=zeros(N,M);  % Kh?i t?o ma tr?n Block 8x8.



k=10;





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
save ('Compressed.txt','I_runcode');
%% PART-2: DECOMPRESSION TECHNIQUE.


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
I_rec= imcrop(I_rec,[0 0 I_dim_old(2) I_dim_old(1)]);










figure,imshow(I_rec,[0,2^n-1]);
figure, imhist(I_rec);
imwrite(I_rec,'Compress.jpg');



% Clearing unused variables from Memory Space.
% clear I_zigzag run;



% Run-length Decoding .

% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.

% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)




% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
