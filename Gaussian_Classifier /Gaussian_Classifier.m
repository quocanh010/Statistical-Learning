S = load('TrainingSamplesDCT_8_new.mat');
Cheetah = S.TrainsampleDCT_FG;
Background = S.TrainsampleDCT_BG;

%Calculating Piror Problability 
Prior_Cheetah = size(Cheetah, 1) / (size(Cheetah,1) + size(Background, 1));
Prior_Background = size(Background, 1) / (size(Cheetah,1) + size(Background, 1));

%Calculate the mean for each feature
mean_Cheetah    = mean(Cheetah, 1);
mean_Background = mean(Background, 1);

%Find covariance matrix for each class
cov_Cheetah = cov(Cheetah);
cov_Background = cov(Background);

%Graph marginal distributions
figure
for i = 1:1:64
    mu_c = mean_Cheetah(i);
    sigma_c = sqrt(cov_Cheetah(i,i));
    mu_b = mean_Background(i);
    sigma_b = sqrt(cov_Background(i,i));
    x_c = linspace(min(Cheetah(:,i)), max(Cheetah(:,i)), 250);
    x_b = linspace(min(Background(:,i)), max(Background(:,i)), 1053);
    
    if (mod((i-1), 8) == 0)
        figure
        count = 1;
    end
    subplot(2, 4, count)
    plot(x_c, pdf('Normal', x_c, mu_c, sigma_c))
    hold on
    plot(x_b, pdf('Normal', x_b, mu_b, sigma_b))
    title(['Feature: ', num2str(i)])
    hold off
    count = count + 1;
    
end
% Draw histogram of the index

% Calculating P_Y_g_X
%img = imread('cheetah.bmp')