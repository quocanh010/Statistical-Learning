S = load('TrainingSamplesDCT_8.mat')
Cheetah = S.TrainsampleDCT_FG
Background = S.TrainsampleDCT_BG

%Calculating Piror Problability 
Prior_Cheetah = size(Cheetah, 1) / (size(Cheetah,1) + size(Background, 1))
Prior_Backgtound = size(Background, 1) / (size(Cheetah,1) + size(Background, 1))

%Find the second largest DCT index
[Cheetah_2nd, I_Cheetah] = max(abs(Cheetah(:,2:end)), [], 2)
[Background_2nd, I_Backgrough] = max(abs(Background(:,2:end)), [], 2)
%Add back the first index
I_Cheetah = I_Cheetah + 1
I_Backgrough = I_Backgrough + 1
% Draw histogram of the index
P_X_g_Y_Cheetah = histogram(I_Cheetah, 'Normalization','probability')
P_X_g_Y_Background = histogram(I_Backgrough , 'Normalization','probability')


% Calculating P_Y_g_X
img = imread('cheetah.bmp')
%compute dct
