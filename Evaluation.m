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
    
    % Uncomment the function you want to evaluate
    %[reconstructed_image, sizeCompressedBits] = jpeg_compression_standard_implementation(original_image, Q);
    [reconstructed_image, sizeCompressedBits] = jpeg_compression_one_encoding(original_image, Q);
    %[reconstructed_image, sizeCompressedBits] = jpeg_compression_seperate_encoding(original_image, Q);

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