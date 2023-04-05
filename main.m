clc
close all
clear

%% input parameters

h_pixel = 50; % height of final pixelated image
w_pixel = 30; % width of final pixelated image
K = 5; % amount of different colors in the final image
source_filename = "obamna.jpg"; % source file name

%% open image

input_imag = imread(source_filename);
input_imag = rgb2lab(input_imag);
figure
imshow(lab2rgb(input_imag))
title('OBAMNA')

%% ---------------- algoritm ---------------------------------

%% initialize superpixels, palette and temperature

N = h_pixel*w_pixel; %set of output pixels
M = prod(size(input_imag(:,:,1))); %get total amount of pixels of input image
input_superpixels = zeros(size(input_imag(:,:,1))); %create array that maps pixels to superpixels
superpixel_colors = zeros(3, N); % give each superpixel a color
temperature = zeros(3,1);
temperature_c = zeros(3,1);
[h_input, w_input] = size(input_superpixels(:,:,1)); % get width and height of input image

%assigning input pixels to superpixels
superpixel_height = round(h_input/h_pixel); %round if wanted image size is not N times input image size
superpixel_width = round(w_input/w_pixel);

counter = 1;
for i =1:superpixel_height:h_input-superpixel_height %assign each input pixel to a superpixel
    for j = 1:superpixel_width:w_input-superpixel_width
        input_superpixels(i:i+superpixel_height, j:j+superpixel_width) = counter;
        counter = counter + 1;
    end
end

L_initial_color = mean(input_imag(:,:,1), 'all'); %get mean L color of input image
a_initial_color = mean(input_imag(:,:,2), 'all'); %get mean L color of input image
b_initial_color = mean(input_imag(:,:,3), 'all'); %get mean L color of input image

superpixel_colors= [repelem(L_initial_color, N); repelem(a_initial_color, N); repelem(b_initial_color, N)]; % give each superpixel the same color

%get temeperature_c
C_L = max(cov(input_imag(:,:,1)), [], 'all'); %get max eigenvalue of input image
C_a = max(cov(input_imag(:,:,2)), [], 'all');
C_b = max(cov(input_imag(:,:,3)), [], 'all');

temperature_c = 2*[C_L C_a C_b]; % set critical temperature
temperature = 1.1*temperature_c; %initial value set for temperature

%% loop for T > Tf

    %% refine superpixels with modified SLIC

    %% associate superpixels to colors in the palette

    %% refine colors in the palette

    %% check convergence of palette

%% post process