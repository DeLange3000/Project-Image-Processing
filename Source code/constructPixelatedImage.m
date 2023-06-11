function [pixelated_image] = constructPixelatedImage(superpixel_palette_colors, palette, h_pixel, w_pixel)

    pixelated_image = zeros(h_pixel, w_pixel, 3);

    for i = 1:h_pixel
        for j = 1:w_pixel
            pixelated_image(i,j, :) = palette(:, superpixel_palette_colors((i-1)*w_pixel + j));
        end
    end
end