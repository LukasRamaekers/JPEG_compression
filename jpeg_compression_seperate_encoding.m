function [Inverse_DCT_rgb_image, compression_size] = jpeg_compression_seperate_encoding(rgb_image, quality)

%==============DCT and quantization==============%

% Convert RGB to YCbCr
ycbcr_image = rgb2ycbcr(rgb_image);
original_size = size(ycbcr_image);

% Calculate the padding needed for both dimensions
pad_rows = 8 - mod(size(rgb_image, 1), 8);
pad_cols = 8 - mod(size(rgb_image, 2), 8);

% Pad the image to the nearest multiple of 8
padded_ycbcr_image = padarray(ycbcr_image, [pad_rows, pad_cols], 0, 'post');

% Extract individual color channels
Y_channel = padded_ycbcr_image(:,:,1);
Cb_channel = padded_ycbcr_image(:,:,2);
Cr_channel = padded_ycbcr_image(:,:,3);

% Apply DCT to each channel
DCT_Y = blockproc(Y_channel, [8 8], @(block) dct2(block.data));
DCT_Cb = blockproc(Cb_channel, [8 8], @(block) dct2(block.data));
DCT_Cr = blockproc(Cr_channel, [8 8], @(block) dct2(block.data));

%Define the quantization matrices
[QY, QC] = generate_quantization_matrix(quality);

%Perform quantization on each channel
Quantized_DCT_Y = blockproc(DCT_Y, [8 8], @(block) round(block.data ./ QY));
Quantized_DCT_Cb = blockproc(DCT_Cb, [8 8], @(block) round(block.data ./ QC));
Quantized_DCT_Cr = blockproc(DCT_Cr, [8 8], @(block) round(block.data ./ QC));

%==============Huffman encoding==============%

%Encode Y, Cb and C vectors
[y_hcode, y_dict, size_layer_Y] = huffman_encoding(Quantized_DCT_Y);
[Cb_hcode, Cb_dict, size_layer_Cb] = huffman_encoding(Quantized_DCT_Cb);
[Cr_hcode, Cr_dict, size_layer_Cr] = huffman_encoding(Quantized_DCT_Cr);

%==============Huffman decoding==============%

% Huffman Decoding
[decoded_quantized_DCT_Y] = huffman_decoding(y_hcode, y_dict, size_layer_Y(1), size_layer_Y(2));
[decoded_quantized_DCT_Cb] = huffman_decoding(Cb_hcode, Cb_dict, size_layer_Cb(1), size_layer_Cb(2));
[decoded_quantized_DCT_Cr] = huffman_decoding(Cr_hcode, Cr_dict, size_layer_Cr(1), size_layer_Cr(2));

%==============inverseDCT and dequantization==============%

% Dequantization for each channel
Inverse_DCT_Y = blockproc(decoded_quantized_DCT_Y, [8 8], @(block) (block.data .* QY));
Inverse_DCT_Cb = blockproc(decoded_quantized_DCT_Cb, [8 8], @(block) (block.data .* QC));
Inverse_DCT_Cr = blockproc(decoded_quantized_DCT_Cr, [8 8], @(block) (block.data .* QC));

% Apply inverse DCT
Inverse_DCT_Y = blockproc(Inverse_DCT_Y, [8 8], @(block) idct2(block.data));
Inverse_DCT_Cb = blockproc(Inverse_DCT_Cb, [8 8], @(block) idct2(block.data));
Inverse_DCT_Cr = blockproc(Inverse_DCT_Cr, [8 8], @(block) idct2(block.data));

% Remove padding from the inverse DCT-transformed channels (may need to be altered since originalm size is usually not known)
Inverse_DCT_Y = Inverse_DCT_Y(1:original_size(1), 1:original_size(2));
Inverse_DCT_Cb = Inverse_DCT_Cb(1:original_size(1), 1:original_size(2));
Inverse_DCT_Cr = Inverse_DCT_Cr(1:original_size(1), 1:original_size(2));

% Combine the inverse DCT-transformed channels back into YCbCr image
Inverse_DCT_ycbcr_image = cat(3, Inverse_DCT_Y, Inverse_DCT_Cb, Inverse_DCT_Cr);

% Convert back to RGB for visualization
Inverse_DCT_rgb_image = ycbcr2rgb(uint8(Inverse_DCT_ycbcr_image));

%{
% Display the original and reconstructed images
figure;
subplot(1, 2, 1), imshow(rgb_image), title('Original Color Image');
subplot(1, 2, 2), imshow(Inverse_DCT_rgb_image), title('Reconstructed Color Image');
%}

% Display evaluation metrics
compression_size = numel(y_hcode) + numel(Cb_hcode) + numel(Cr_hcode);
end

%==============Quantization quality==============%

function [QY, QC] = generate_quantization_matrix(quality)
    % Standard JPEG quantization matrix for luminance (Y) channel
    QY = [
    16  11  10  16  24  40  51  61;
    12  12  14  19  26  58  60  55;
    14  13  16  24  40  57  69  56;
    14  17  22  29  51  87  80  62;
    18  22  37  56  68 109 103  77;
    24  35  55  64  81 104 113  92;
    49  64  78  87 103 121 120 101;
    72  92  95  98 112 100 103  99
    ];
    
    QY = scale_quantization_table(QY, quality);
    
    % Standard JPEG quantization matrix for chrominance (Cb and Cr) channels
    QC = [
    17  18  24  47  99  99  99  99;
    18  21  26  66  99  99  99  99;
    24  26  56  99  99  99  99  99;
    47  66  99  99  99  99  99  99;
    99  99  99  99  99  99  99  99;
    99  99  99  99  99  99  99  99;
    99  99  99  99  99  99  99  99;
    99  99  99  99  99  99  99  99
    ];

    QC = scale_quantization_table(QC, quality);
end

function scaled_table = scale_quantization_table(base_table, quality)
    % Ensure quality is within the valid range
    quality = max(1, min(100, quality));

    % Compute the scaling factor S based on the quality factor Q
    S = (quality < 50) .* (5000 ./ quality) + (quality >= 50) .* (200 - 2 * quality);

    % Initialize the scaled table
    scaled_table = zeros(size(base_table));

    % Calculate the scaled values using the formula T_s[i] = floor((S * T_b[i] + 50) / 100)
    for i = 1:numel(base_table)
        scaled_value = floor((S * base_table(i) + 50) / 100);

        % Ensure that the scaled value is at least 1
        scaled_value = max(1, scaled_value);

        % Set the scaled value in the table
        scaled_table(i) = scaled_value;
    end
end


%==============Huffman algorithm==============%

function [hcode, dict, size_layer] = huffman_encoding(channel)
    % Combine Y, Cb, and Cr into a single vector
    vector = zigzag(channel);
   
    size_layer = size(channel);

    % Size of the combined vector
    [m, n] = size(vector);
    Totalcount = m * n;
    
    % Variables to find the probability
    cnt = 1;
    sigma = 0;
    
    % Computing the cumulative probability.
    for i = min(vector):max(vector)
        k = vector == i;
        count(cnt) = sum(k(:));
        
        % Pro array is having the probabilities
        pro(cnt) = count(cnt) / Totalcount;
        sigma = sigma + pro(cnt);
        cnt = cnt + 1;
    end
    
    % Normalize the probability vector
    pro = pro / sum(pro);

    % Symbols for an image
    symbols = min(vector):max(vector);

    % Use the probability vector for Huffman dictionary creation
    dict = huffmandict(symbols, pro);
    
    % Huffman Encoding
    hcode = huffmanenco(vector, dict);
end


function [decoded_channel] = huffman_decoding(hcode, dict, m, n)
    % Huffman Decoding
    decoded_vector = huffmandeco(hcode, dict);

    % Return the channel back to its original size
    decoded_channel = izigzag(decoded_vector, m, n);

end
