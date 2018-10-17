S = load('TrainingSamplesDCT_8.mat')
Cheetah = S.TrainsampleDCT_FG
Background = S.TrainsampleDCT_BG

%Calculating Piror Problability 
Prior_Cheetah = size(Cheetah, 1) / (size(Cheetah,1) + size(Background, 1))
Prior_Background = size(Background, 1) / (size(Cheetah,1) + size(Background, 1))

%Find the second largest DCT index
[~, I_Cheetah] = max(abs(Cheetah(:,2:end)), [], 2)
[~, I_Backgrough] = max(abs(Background(:,2:end)), [], 2)
%Add back the first index
I_Cheetah = I_Cheetah + 1
I_Backgrough = I_Backgrough + 1
% Draw histogram of the index
P_X_g_Y_Cheetah = histogram(I_Cheetah,'NumBins', 63, 'BinLimits', [2 64], 'Normalization','probability')
C_Cheetah =P_X_g_Y_Cheetah.Values
C_Cheetah = [0, C_Cheetah]
P_X_g_Y_Background = histogram(I_Backgrough,'NumBins', 63, 'BinLimits', [2 64], 'Normalization','probability')
C_Background = P_X_g_Y_Background.Values
C_Background = [0, C_Background]
% Calculating P_Y_g_X
img = imread('cheetah.bmp')
% Normalizing Image
img = im2double(img)

%Vectorize 8-8 block with zig-zag pattern
fileID = fopen('Zig-Zag Pattern.txt')
formatSpec = '%i'
index_v = fscanf(fileID,formatSpec)
index_v = int16(index_v)
index_v = index_v + 1
%Padding to teh image
B = padarray(img, [7 7], 'symmetric','post')
%compute ans sliding window
%  for i = 1:size(img, 1) 
%      for j = 1:size(img, 2) 
%          result_v(count) = zig_zag_v(B(i:i+7, j:j+7), index_v, C_Cheetah, C_Background, Prior_Cheetah, Prior_Background)
%          count = count + 1
%      end
%  end
result = nlfilter(  B , [8 8], @zig_zag_v)
Cheetah_m = C_Cheetah(result) .* Prior_Cheetah
Background_m = C_Background(result) .* Prior_Background
n_img = Background_m <= Cheetah_m
n_img = n_img(1:255, 1:270)
imagesc(n_img)
colormap(gray(256))
axis equal

masked_cheetah = imread('cheetah_mask.bmp')
