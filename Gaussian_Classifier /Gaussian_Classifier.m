S = load('TrainingSamplesDCT_8_new.mat');
Cheetah = S.TrainsampleDCT_FG;
Background = S.TrainsampleDCT_BG;

%Calculating Piror Problability 
Prior_Cheetah = size(Cheetah, 1) / (size(Cheetah,1) + size(Background, 1));
Prior_Background = size(Background, 1) / (size(Cheetah,1) + size(Background, 1));

%Find the second largest DCT index
[~, I_Cheetah] = max(abs(Cheetah(:,2:end)), [], 2);
[~, I_Backgrough] = max(abs(Background(:,2:end)), [], 2);
%Add back the first index
I_Cheetah = I_Cheetah + 1;
I_Backgrough = I_Backgrough + 1;
% Draw histogram of the index

% Calculating P_Y_g_X
img = imread('cheetah.bmp')