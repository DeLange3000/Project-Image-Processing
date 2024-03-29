clc
close all
clear

%% input parameters

h_pixel = 30; % height of final pixelated im age
w_pixel = 40; % width of final pixelated image
K = 8; % amount of different colors in the final image
perturbation_factor = 50;
addpath('images');
% source_filename = "obamna.jpg"; % source file name
% source_filename = "shrek.jpg"; % source file name
source_filename = "dolphin.jpg"; % source file name


%% open image

input_imag = imread(source_filename);
input_imag = rgb2lab(input_imag);
input_imag = imresize(input_imag, 0.5);
figure
imshow(lab2rgb(input_imag))
title('Original Image')
hold on

%% ---------------- algorithm ---------------------------------

%% initialize superpixels, palette and temperature

N = h_pixel*w_pixel; %set of output pixels, number of superpixels
M = numel(input_imag(:,:,1)); %get total amount of pixels of input image
input_superpixels = zeros(size(input_imag(:,:,1))); %create array that maps pixels to superpixels. Each pixel 'index' will have int value in
                                                    %this array that corresponds to the superpixel they belong to
superpixel_colors = zeros(3, N); % give each superpixel a color
superpixel_palette_colors = round(rand(N,1))+1; %for assigning random palette color to superpixels (points to palette index)
temperature = zeros(3,1);
temperature_c = zeros(3,1);
temperature_f = 1;
output_imag = zeros(h_pixel, w_pixel, 3);
[h_input, w_input] = size(input_superpixels(:,:,1)); % get width and height of input image

% palette = zeros(3,K);
b_and_a_sat = 1.1;


%assigning input pixels to superpixels
superpixel_height = h_input/h_pixel; %round if wanted image size is not N times input image size
superpixel_width = w_input/w_pixel;
superpixel_center = zeros(2, N);

counter = 1;
for i = 1:superpixel_height:h_input %assign each input pixel to a superpixel
    for j = 1:superpixel_width:w_input
        jr = round(j);
        ir = round(i);

        input_superpixels(ir:ir+superpixel_height + 1, jr:jr+superpixel_width + 1) = counter;
        counter = counter + 1;
    end
end

input_superpixels = input_superpixels(1:h_input, 1:w_input); % get correct array size

% get superpixel colors
superpixel_colors = zeros(N, 3);
for i = 1:N
    [x, y] = find(input_superpixels == i);
    superpixel_colors(i, :) = mean(input_imag(x,y, :), [1, 2]);
end

%get temperature_c (idk if correct)
coefficients_L = pca(input_imag(:,:,1));
PC_L = coefficients_L(:, 1);
C_L = var(input_imag(:,:,1) * PC_L);
coefficients_a = pca(input_imag(:,:,2));
PC_a = coefficients_a(:, 1);
C_a = var(input_imag(:,:,2) * PC_a);
coefficients_b = pca(input_imag(:,:,3));
PC_b = coefficients_b(:, 1);
C_b = var(input_imag(:,:,3) * PC_b);

% imagesc(input_superpixels)
initial_color = squeeze((mean(input_imag, [1, 2]))); %get mean L color of input image
p_ck = [0.5 0.5]; %initiate P(ck)
p_ps = 1/N;
palette = [initial_color, initial_color];
e_palette = 5; % maximum change for new color
e_cluster = 1; %minimum change for palette color split
temperature_lowering_factor = 0.7; %0.7;

%subclusters
subclusters = [1; 2];

% C_L = max(eig(cov(input_imag(:,:,1)))); %get max eigenvalue of input image
% C_a = max(eig(cov(input_imag(:,:,2))));
% C_b = max(eig(cov(input_imag(:,:,3))));

temperature_c = 2*norm([C_L C_a C_b]); % set critical temperature
temperature = 1.1*temperature_c; %initial value set for temperature

% get center of each superpixel (assume they are rather square)
for i = 1:N
    [x,y] = find(input_superpixels == i); %find elements equal to i
    x_superpixel_centre = (max(x) + min(x))/2;
    y_superpixel_centre = (max(y) + min(y))/2;
    superpixel_center(:,i) = [y_superpixel_centre, x_superpixel_centre];
end
%plot(superpixel_center(1,:),superpixel_center(2,:), "*r")


%% loop for T > Tf

loop_index = 0;
while temperature > temperature_f
    disp('temperature: '+string(temperature))
    loop_index = loop_index + 1;
    disp('loop:' + string(loop_index))

    %% refine superpixels with modified SLIC (simple linear iterative clustering)

    % euclidean distance = sqrt((x1-x2)^2 + (y1-y2)^2 + ...)
    m = 45;
    total_palette_change = palette;
    if(length(palette(1,:)) ~= N)
        for i = 1:length(palette(1,:))
            palette(:,i) = palette(:,i) + [rand; randn; randn].*perturbation_factor; %permutations of mean
        end
    else
        e_palette = 1;
         palette(:,i) = palette(:,i) + [rand; randn; randn].*perturbation_factor/100;
    end

    %iterate for each input pixel over superpixels that it is next too (got to optimize
    %this)
            minimal_d = ones(h_input, w_input)*1e9;
            
            if (loop_index >1)
                superpixel_colors = msprime_superpixels';
            else
                superpixel_colors = superpixel_colors';
            end
            for aa = 1:N
			    %color
                pc_input = input_imag;
                pc_super = palette(:, superpixel_palette_colors(aa)); %superpixel_colors(:,aa); % should be superpixel_palette_colors(aa)
                a1 =repmat(pc_super(1), h_input, w_input);
                a2 =repmat(pc_super(2), h_input, w_input);
                a3 =repmat(pc_super(3), h_input, w_input);
                pc_super =  cat(3, a1, a2, a3);
                dc = sqrt(sum((pc_input - pc_super).^2, 3));
        
                %position
                pd_input = cat(3,ones(h_input, w_input).*[1:h_input]', ones(h_input, w_input).*[1:w_input]);
                pd_super = cat(3, superpixel_center(2,aa)*ones(h_input, w_input), superpixel_center(1,aa)*ones(h_input, w_input));
                dp = m.*sqrt(N/M).*sqrt((pd_input(:,:,1) - pd_super(:,:,1)).^2 + (pd_input(:,:,2) - pd_super(:,:,2)).^2);
    
			    d_min = dc+dp;
                change_indices = d_min < minimal_d;
                minimal_d(change_indices) = d_min(change_indices);
                input_superpixels(change_indices) = aa;
            end
% 	figure
% 	imagesc(input_superpixels)
% 	hold on
    %% laplacian smoothing (hopefully doesnt cause mega pixel bugs anymore)

    for i = (w_pixel+1):(N-w_pixel+1) % skip edges
        %get neigbouring superpixels
        if (mod(i, w_pixel) == 1 || mod(i, w_pixel) == 0) %skip left and right border of image
            neigboring_superpixels = [];
        else
            neigboring_superpixels = [i - 1, i + 1, i + w_pixel, i - w_pixel];
            mean_position = mean(superpixel_center(:,neigboring_superpixels), 2);
            superpixel_center(:,i) = 0.6*superpixel_center(:,i) + 0.4*mean_position; % move superpixel centers
        end
	    
    end
    
    % smoothing color representatives of superpixels Original SLIC
    % algorithm represents every superpixel using ms as color, this is the
    % mean color of all pixels in superpixel. This can give problems with
    % continuous regions (color gradient). This region will not seem
    % continuous in the pixelated output. To solve this, the values for ms
    % are altered by using a bilateral filter (=non-linear,
    % edge-preserving, and noise-reducing smoothing filter for images).
    % This gives the new colors => ms'
    
    %==========Map every superpixel to output grid==========
    
    %Get mean color of every superpixel
    size_image = size(input_imag);
    size_superpixels = size(input_superpixels);

    ms_superpixels = zeros(N, 3); %ms value for every superpixel
    superpixels_image = zeros(size(input_imag));
    % Iterate over each superpixel
    for i = 1:N
        % Find the indices in input_superpixels that have the value
        % corresponding to superpixel i
        [row_indices, col_indices] = find(input_superpixels == i);
        
        if (not(isempty(row_indices))) % only run for still existing superpixels
            color_values = input_imag(row_indices, col_indices, :); % Get the corresponding color channel values in input_imag
			ms_superpixels(i,:) = mean(color_values, [1 2]);% Calculate the mean of the corresponding color channel values in input_imag
        end
    end
    
    %show image with all pixels set to their superpixel color
    for i=1:h_input
        for j=1:w_input
            superpixels_image(i,j,:) = ms_superpixels(input_superpixels(i,j),:);
        end
    end
%      figure
%      imshow(lab2rgb(superpixels_image))
%      hold on
%      %plot(superpixel_center(1,:),superpixel_center(2,:), "*r")
%      title('superpixels of image')
    
    %==========Convert superpixels to grid==========
    %Create output image based on the superpixels
    %nie zeker of da zo moet, werkt alleen als de volgorde van de
    %superpixels = de volgorde van de output pixels
    output_pixels = ms_superpixels;
    output_pixels = reshape(output_pixels,w_pixel,h_pixel,3);
    output_pixels = permute(output_pixels,[2 1 3]);    
    
    % Bilateral filtering of this image
    filtered_image = imbilatfilt(lab2rgb(output_pixels));
    filtered_image = rgb2lab(filtered_image);
    msprime_superpixels = permute(filtered_image,[2 1 3]);
    msprime_superpixels = reshape(msprime_superpixels,N,3);

%     figure
%     imshow([lab2rgb(output_pixels) lab2rgb(filtered_image)],'InitialMagnification',2000)
%     title('Before and after bilaterial filtering')

    
    %% refine colors in the palette

    %calculate p(ck|ps)

    p_ck_ps = zeros(h_pixel, w_pixel, width(palette));
    for i = 1:h_pixel
        for j = 1:w_pixel
            normalisation_factor = GetNormalisationFactor(p_ck, msprime_superpixels((i-1)*w_pixel + j, :), palette, temperature); 
            for a = 1:width(palette)
                p_ck_ps(i,j,a) = p_ck(a)*exp(-norm(msprime_superpixels((i-1)*w_pixel + j, :) - palette(:,a)')/temperature)/normalisation_factor;
            end
            palette_color_index = find(p_ck_ps(i,j,:) == max(p_ck_ps(i,j,:))); % gets indexes that points to which color in the palette is most probable
            %assign random color to superpixel if 50/50% chance
            if(length(palette_color_index) > 1)
                superpixel_palette_colors((i-1)*w_pixel+j) = palette_color_index(randi(length(palette_color_index))); % assign palette color to superpixel
            else
                 superpixel_palette_colors((i-1)*w_pixel+j) = palette_color_index;
            end

        end
    end

    p_ck = zeros(width(palette),1);
    for a = 1:width(palette)
         p_ck(a) = sum(p_ck_ps(:,:,a).*p_ps, "all");
    end

    %refine colors
    temp_palette = zeros(3,length(palette(1,:)));
    for a = 1:width(palette)
        for i = 1:h_pixel
            for j = 1:w_pixel
                temp_palette(:,a) = temp_palette(:,a) + msprime_superpixels((i-1)*w_pixel + j, :)'.*p_ck_ps(i,j,a).*p_ps/p_ck(a);
            end
        end
    end

    palette = temp_palette;
    disp('palette: '+string(palette))
    %disp('P_ck_ps: '+string(p_ck_ps))


    %% check convergence of palette

    %check T
    total_palette_change = norm(total_palette_change - palette);

    if(total_palette_change < e_palette)
        temperature = temperature*temperature_lowering_factor;
        if(length(palette(1,:)) < K)
                for i = 1:length(palette(1,:)) - 1
                    ck1 = palette(:,subclusters(1,i));
                    ck2 = palette(:,subclusters(2,i));
                    if(norm(ck1 - ck2) > e_cluster)
                        palette = [palette ck1 + 100*[C_L; C_a; C_b]./temperature_c];
                        palette(:, i) = ck2 - 100*[C_L; C_a; C_b]./temperature_c;
                        subclusters = [subclusters [i; length(palette)]];
                        p_ck = [p_ck; p_ck(i)/2];
                        p_ck(i) = p_ck(i)/2;
                        break
                    end

                end
        end

    end
    % show temp image
    pixelated_image = constructPixelatedImage(superpixel_palette_colors, palette, h_pixel, w_pixel);
 
%     figure
%     imshow(lab2rgb(pixelated_image),'InitialMagnification',400)
%     title('Pixelated Image')
end
%% post process

palette(2:3,:) = b_and_a_sat*palette(2:3,:);

%% show image

pixelated_image = constructPixelatedImage(superpixel_palette_colors, palette, h_pixel, w_pixel);
 
figure
imshow(lab2rgb(pixelated_image),'InitialMagnification',2000)
title('Pixelated Image')
