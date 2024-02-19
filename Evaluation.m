% Specify the path to the image
[filename, filepath] = uigetfile({'*.jpg;*.png;*.bmp;*.tif', 'Image files (*.jpg, *.png, *.bmp, *.tif)'; '*.*', 'All files (*.*)'}, 'Select an Image');

% Load the image
original_image = imread(fullfile(filepath, filename));

info = imfinfo(fullfile(filepath, filename));

[Nx, Ny] = size(original_image);

fprintf('\n'); 
fprintf('Format:   %s',info.Format); fprintf('\n'); 
fprintf('Size:     Width = %d,  Height = %d',info.Width,info.Height); fprintf('     [pixels] \n\n');

% Size of the original image in bits
sizeOriginalBits = info.Width * info.Height * 8;

k = 0;
compressionRatios = zeros(1, 98);
Rate = zeros(1, 98);
Distortion = zeros(1, 98);
PSNR = zeros(1, 98);

for Q = 0:1:97
    k = k+1    
    [reconstructed_image, sizeCompressedBits] = jpeg_compression_seperate_encoding(original_image, Q);

    % Calculate compression ratio
    compressionRatios(k) = sizeOriginalBits / sizeCompressedBits;

    Rate(k) = sizeCompressedBits/(Nx*Ny); 
    
    Distortion(k) = immse(double(reconstructed_image),double(original_image));
    PSNR(k) = psnr(double(reconstructed_image),double(original_image),255);
end

figure(1)
subplot(2, 2, 1), imshow(original_image), title('Original Image');

subplot(2, 2, 2), plot(Rate,Distortion); grid on;
title('JPEG rate distortion curve');
xlabel('Rate  [bpp]');
ylabel('MSE');

subplot(2, 2, 3), plot(Rate,PSNR); grid on;
title('JPEG rate distortion curve');
xlabel('Rate  [bpp]');
ylabel('PSNR  [dB]');

subplot(2, 2, 4), plot(0:97, compressionRatios, '-o');
grid on;
title('Compression Ratio vs. JPEG Quality');
xlabel('JPEG Quality');
ylabel('Compression Ratio');

%{
%==============Evaluation==============%
function evaluate_compression(original_image, reconstructed_image, compressed_size)
    % Calculate compression ratio
    original_size = numel(original_image);
    compression_ratio = original_size / compressed_size;
    
    % Calculate rate
    rate = compressed_size / original_size;
    
    % Calculate distortion (Mean Squared Error)
    mse = immse(double(reconstructed_image),double(original_image));
    
    % Calculate PSNR
    psnr_value = psnr(double(reconstructed_image),double(original_image),255);
    
    % Display the results
    fprintf('Compression Ratio: %.4f\n', compression_ratio);
    fprintf('Rate: %.4f\n', rate);
    fprintf('Distortion (MSE): %.4f\n', mse);
    fprintf('PSNR: %.4f dB\n', psnr_value);
end
%}