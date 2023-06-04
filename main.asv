clc
close all
clear

%% input parameters

h_pixel = 10; % height of final pixelated image
w_pixel = 10; % width of final pixelated image
K = 5; % amount of different colors in the final image
loop_range = 100; % how many iterations for image refinement
addpath('images');
% source_filename = "obamna.jpg"; % source file name
% source_filename = "shrek.jpg"; % source file name
source_filename = "images/Two-colours-wallpaper1.jpg"; % source file name


%% open image

input_imag = imread(source_filename);
input_imag = rgb2lab(input_imag);
input_imag = imresize(input_imag, 0.5);
figure
imshow(lab2rgb(input_imag))
title('OBAMNA')
hold on

%% ---------------- algorithm ---------------------------------

%% initialize superpixels, palette and temperature

N = h_pixel*w_pixel; %set of output pixels, number of superpixels
M = numel(input_imag(:,:,1)); %get total amount of pixels of input image
input_superpixels = zeros(size(input_imag(:,:,1))); %create array that maps pixels to superpixels. Each pixel 'index' will have int value in
                                                    %this array that corresponds to the superpixel they belong to
superpixel_colors = zeros(3, N); % give each superpixel a color
superpixel_palette_colors = ones(N, 1); %for assigning palette colors to superpixels (points to palette index)
temperature = zeros(3,1);
temperature_c = zeros(3,1);
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
        input_superpixels(ir:ir+superpixel_height, jr:jr+superpixel_width) = counter;
        counter = counter + 1;
    end
end
% imagesc(input_superpixels)

initial_color = squeeze((mean(input_imag, [1, 2]))); %get mean L color of input image
p_ck = [1]; %initiate P(ck)
p_ps = 1/N;
palette = [initial_color];
e_palette = 300; % maximum change for new color
e_cluster = 10;
temperature_lowering_factor = 0.7;

%subclusters??????????????
sub_disturbance = 30;
subclusters = [1; 2];

%get temperature_c (idk if correct)
C_L = max(eig(cov(input_imag(:,:,1)))); %get max eigenvalue of input image
C_a = max(eig(cov(input_imag(:,:,2))));
C_b = max(eig(cov(input_imag(:,:,3))));

temperature_c = 2*norm([C_L C_a C_b]); % set critical temperature
temperature = 1.1*temperature_c; %initial value set for temperature

% get center of each superpixel (assume they are rather square)
for i = 1:N
    minx = w_input;
    maxx = 0;
    miny = h_input;
    maxy = 0;
    [x,y] = find(input_superpixels == i); %find elements equal to i
    x_superpixel_centre = (max(x) + min(x))/2;
    y_superpixel_centre = (max(y) + min(y))/2;
    superpixel_center(:,i) = [y_superpixel_centre, x_superpixel_centre];
end
plot(superpixel_center(1,:),superpixel_center(2,:), "*r")

%% loop for T > Tf

for loop_index = 1:loop_range
    disp(loop_index)

    %% refine superpixels with modified SLIC (simple linear iterative clustering)

    % euclidean distance = sqrt((x1-x2)^2 + (y1-y2)^2 + ...)
    m = 45;
    total_palette_change = palette;



    %iterate for each input pixel over superpixels that it is next too (got to optimize
    %this)
    for i = 1:h_input
        for j = 1:w_input
%             if(mod(i*j, 10000) == 0)
%                 disp("looping... "+num2str(i*j)+" of "+num2str(h_input*w_input))
%             end
            
            superpixel = input_superpixels(i,j);
            new_superpixel = superpixel;
         
            %get neigboring superpixels
            if (mod(superpixel, w_pixel) == 1) %left border of image
                neigboring_superpixels = [superpixel, superpixel+1];
            elseif (mod(superpixel, w_pixel) == 0) %right border of image
                neigboring_superpixels = [superpixel - 1, superpixel];
            else %center of image
                neigboring_superpixels = [superpixel - 1, superpixel, superpixel+1];
            end
			
            if (superpixel <= w_pixel) %top border of image
                neigboring_superpixels = [neigboring_superpixels, neigboring_superpixels + w_pixel];
			elseif superpixel > N - w_pixel %bottom border of image
                neigboring_superpixels = [neigboring_superpixels - w_pixel, neigboring_superpixels];
            else %center of image
                neigboring_superpixels = [neigboring_superpixels - w_pixel, neigboring_superpixels, neigboring_superpixels + w_pixel];
			end


			%disp(neigboring_superpixels)
			
			% dmin has to first be (i,j)'s current superpixel, otherwise it
			% will by default be moved to the leftmost neighboring
			% superpixel.
			%color
            pc_input = input_imag(i,j,:);
            pc_super = palette(:, superpixel_palette_colors(superpixel));
            dc = sqrt(sum((pc_input - pc_super).^2, 'all'));
    
            %position
            pd_input = [i,j];
            pd_super = [x_superpixel_centre y_superpixel_centre];
            dp = m*sqrt(N/M)*sqrt(sum((pd_input - pd_super).^2));

			d_min = dc+dp;

            % loop over neighboring superpixels
            for a = neigboring_superpixels

                %color
                pc_input = input_imag(i,j,:);
                pc_super = palette(:, superpixel_palette_colors(a));
                dc = sqrt(sum((pc_input - pc_super).^2, 'all'));
    
                %position
                pd_input = [i,j];
                pd_super = [x_superpixel_centre y_superpixel_centre];
                dp = m*sqrt(N/M)*sqrt(sum((pd_input - pd_super).^2));
				

                if(dc + dp < d_min) %check euclidean distance
                    d_min = dc + dp;
                    new_superpixel = a;
                end
			end
            input_superpixels(i,j) = new_superpixel; %assign pixel to new superpixel
			
        end
	end
% 	figure
% 	imagesc(input_superpixels)
% 	hold on
    %% laplacian smoothing

    % get center of each superpixel (assume they are rather square) (copied
    % from code above)
    for i = 1:N
        minx = w_input;
        maxx = 0;
        miny = h_input;
        maxy = 0;
        [x,y] = find(input_superpixels == i); %find elements equal to i
        x_superpixel_centre = (max(x) + min(x))/2;
        y_superpixel_centre = (max(y) + min(y))/2;
		if size([y_superpixel_centre, x_superpixel_centre]) ~= [0 2]
			superpixel_center(:,i) = [y_superpixel_centre, x_superpixel_centre];
		end
	end
% 	plot(superpixel_center(1,:),superpixel_center(2,:), "*r")

    % for i = w_pixel+1:N-w_pixel % skip edges
    %     %get neigbouring superpixels
    %     if (mod(i, w_pixel) == 1 || mod(i, w_pixel) == 0) %skip left and right border of image
    %         neigboring_superpixels = [];
    %     else
    %         neigboring_superpixels = [i - 1, i + 1, i + w_pixel, i - w_pixel];
    %     end
	% 
    %     net_distance = [0; 0]; %x and y
    %     for a = neigboring_superpixels
    %         superpixel_center(:,i) - superpixel_center(:,a);
    %         net_distance = net_distance + superpixel_center(:,i) - superpixel_center(:,a);
    %     end
    %     superpixel_center(:,i) = superpixel_center(:,i) + 0.4*net_distance; % move superpixel centers
    % end
    
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

    % compare sizes, input_superpixels was one row too much, fix somewhere else
    if size_image(1) < size_superpixels(1)
        input_superpixels = input_superpixels(1:size_image(1), :);
    end

    if size_image(2) < size_superpixels(2)
        input_superpixels = input_superpixels(:, 1:size_image(2));
    end

    ms_superpixels = zeros(N, 3); %ms value for every superpixel
    superpixels_image = zeros(size(input_imag));

    % Iterate over each superpixel
    for i = 1:N
        % Find the indices in input_superpixels that have the value
        % corresponding to superpixel i
        [row_indices, col_indices] = find(input_superpixels == i);
        
        if (not(isempty(row_indices))) % only run for still existing superpixels
            % Loop over each color channel
            for k = 1:3
                % Get the corresponding color channel values in input_imag
                color_values = input_imag(sub2ind(size(input_imag), row_indices, col_indices, k*ones(size(row_indices))));
                % Calculate the mean of the corresponding color channel values in input_imag
				ms_superpixels(i,k) = mean(color_values);
            end
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

    temp_p_ck = zeros(size(p_ck));
    
    temp_ck = zeros(size(palette));
    for i = 1:h_pixel
        for j = 1:w_pixel
            max_p_ck_ps = 0;
            palette_index = 0;
            p_ck_ps = zeros(1, length(palette(1,:)));
            for a = 1:length(p_ck_ps)
                p_ck_ps(a) = p_ck(a)*exp(-norm(squeeze(filtered_image(i,j,:)) - palette(:,a))/temperature); %needs to be normalized?
                if(p_ck_ps(a) >= max_p_ck_ps)
                    max_p_ck_ps = p_ck_ps;
                    palette_index = a;
                end
            end
            if(sum(p_ck.*exp(-norm(msprime_superpixels((i-1)*w_pixel+j, :)' - palette)./temperature)) ~= 0)
                p_ck_ps = p_ck_ps/sum(p_ck.*exp(-norm(msprime_superpixels((i-1)*w_pixel+j, :)' - palette)./temperature));
            else
                  p_ck_ps = 0;
            end
            for a = 1:length(temp_p_ck)
                temp_p_ck(a) = temp_p_ck(a) + p_ck_ps(a)*p_ps;
                temp_ck(:, a) = temp_ck(:, a) + msprime_superpixels((i-1)*w_pixel+j, :)'*p_ck_ps(a)*p_ps;
            end
            superpixel_palette_colors((i-1)*w_pixel+j) = palette_index;
        end
    end
    p_ck = temp_p_ck;

    %refine colors

    palette = temp_ck./p_ck;


    %% check convergence of palette

    %check T
    if (loop_index ~= loop_range)
    total_palette_change = norm(total_palette_change - palette);

    if(total_palette_change < e_palette)
        temperature = temperature*temperature_lowering_factor;
        if(length(palette(1,:)) < K)
            if(length(palette(1,:)) == 1)
                 palette = [palette*1/3 2/3*palette(:,subclusters(1,1))];
                 p_ck = [p_ck/2 p_ck/2];
            else
                for i = 1:length(palette(1,:)) - 1
                    ck1 = palette(:,subclusters(1,i));
                    ck2 = palette(:,subclusters(2,i));
    
                    if(norm(ck1 - ck2) > e_cluster)
                        palette = [palette ck2 - ck1/2];
                        palette(:, i) = ck2 + ck1/2;
                        subclusters = [subclusters [i; length(palette)]];
                        p_ck = [p_ck p_ck(i)/2];
                        p_ck(i) = p_ck(i)/2;
                    end
                end
            end
            
        end

    else
        palette = palette + 30*randn(size(palette));
    end

    %disturb ck1 and ck2

    % show temp image
    pixelated_image = constructPixelatedImage(superpixel_palette_colors, palette, h_pixel, w_pixel);
 
    figure
    imshow(lab2rgb(pixelated_image),'InitialMagnification',2000)
    title('Pixelated Image')
    end

end
%% post process

palette(2:3,:) = b_and_a_sat*palette(2:3,:);

%% show image

pixelated_image = constructPixelatedImage(superpixel_palette_colors, palette, h_pixel, w_pixel);
 
figure
imshow(lab2rgb(pixelated_image),'InitialMagnification',2000)
title('Pixelated Image')
