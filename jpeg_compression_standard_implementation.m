function [reconstructed_image, sizeCompressedBits] = jpeg_compression_standard_implementation(rgb_image, Q)

imwrite(rgb_image,'X.jpeg','jpg','Quality',Q);
reconstructed_image = imread('X.jpeg');
    
% Size of the compressed image in bits
infoCompressed = imfinfo('X.jpeg');
sizeCompressedBits = infoCompressed.FileSize * 8;