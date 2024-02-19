function [Inverse_DCT_rgb_image, compression_size] = jpeg_compression_one_encoding(rgb_image, quality)

    %==============DCT and quantization==============%
    
    % Convert RGB to YCbCr
    ycbcr_image = rgb2ycbcr(rgb_image);
    original_size = size(ycbcr_image);
    
    % Calculate the padding needed for both dimensions
    pad_rows = 8 - mod(size(ycbcr_image, 1), 8);
    pad_cols = 8 - mod(size(ycbcr_image, 2), 8);
    
    % Pad the image to the nearest multiple of 8
    padded_ycbcr_image = padarray(ycbcr_image, [pad_rows, pad_cols], 0, 'post');
    
    % Extract individual color channels
    Y_channel = padded_ycbcr_image(:,:,1);
    Cb_channel = padded_ycbcr_image(:,:,2);
    Cr_channel = padded_ycbcr_image(:,:,3);
    
    % Define DCT matrix
    dct_matrix = dctmtx(8);
    
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
    [ycbcr_hcode, dict, size_layer] = huffman_encoding_ycbcr(Quantized_DCT_Y, Quantized_DCT_Cb, Quantized_DCT_Cr);
    
    %==============Huffman decoding==============%
    
    % Huffman Decoding
    [decoded_quantized_DCT_Y, decoded_quantized_DCT_Cb, decoded_quantized_DCT_Cr] = huffman_decoding_ycbcr(ycbcr_hcode, dict, size_layer(1), size_layer(2));
    
    %==============inverseDCT and dequantization==============%
    
    % Inverse DCT and dequantization for each channel
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
    
    % Clip values to the valid range [0, 255]
    Inverse_DCT_ycbcr_image = max(0, min(255, Inverse_DCT_ycbcr_image));
    
    % Convert back to RGB for visualization
    Inverse_DCT_rgb_image = ycbcr2rgb(uint8(Inverse_DCT_ycbcr_image));
    
    %{
    % Display the original and reconstructed images
    figure;
    subplot(1, 2, 1), imshow(rgb_image), title('Original Color Image');
    subplot(1, 2, 2), imshow(Inverse_DCT_rgb_image), title('Reconstructed Color Image');
    %}

    % Get the compression size for evaluation
    compression_size = numel(ycbcr_hcode);

end

%==============Quantization quality==============%

function [QY, QC] = generate_quantization_matrix(quality)
    % Quality should be between 1 and 100
    quality = max(1, min(100, quality));
    
    % Adjustment factor
    alpha = 50 / quality;

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
    
    % Adjust the quantization matrix based on the quality parameter
    QY = round((alpha * QY + 50) / 100);
    
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

    % Adjust the quantization matrix for chrominance
    QC = round((alpha * QC + 50) / 100);
end

%==============Huffman algorithm==============%

function [hcode, dict, size_layer] = huffman_encoding_ycbcr(Y, Cb, Cr)
    % Combine Y, Cb, and Cr into a single vector
    zigzag_Y = zigzag(Y);
    zigzag_Cb = zigzag(Cb);
    zigzag_Cr = zigzag(Cr);
    combined_vector = [zigzag_Y(:); zigzag_Cb(:); zigzag_Cr(:)];
    
    size_layer = size(Y);

    % Size of the combined vector
    [m, n] = size(combined_vector);
    Totalcount = m * n;
    
    % Variables to find the probability
    cnt = 1;
    sigma = 0;
    
    % Computing the cumulative probability.
    for i = min(combined_vector):max(combined_vector)
        k = combined_vector == i;
        count(cnt) = sum(k(:));
        
        % Pro array is having the probabilities
        pro(cnt) = count(cnt) / Totalcount;
        sigma = sigma + pro(cnt);
        cumpro(cnt) = sigma;
        cnt = cnt + 1;
    end
    
    % Normalize the probability vector
    pro = pro / sum(pro);

    % Symbols for an image
    symbols = min(combined_vector):max(combined_vector);

    % Use the probability vector for Huffman dictionary creation
    dict = huffmandict(symbols, pro);
    
    % Huffman Encoding
    hcode = huffmanenco(combined_vector, dict);
end

function [Y, Cb, Cr] = huffman_decoding_ycbcr(ycbcr_hcode, dict, m, n)
    % Huffman Decoding
    combined_vector_dec = huffmandeco(ycbcr_hcode, dict);

    % Split the combined vector back into Y, Cb, and Cr
    Y_size = m * n;
    Cb_size = m * n;
    Cr_size = m * n;
    
    % Extract the zigzagged components from the combined vector
    zigzag_Y_dec = combined_vector_dec(1:Y_size);
    zigzag_Cb_dec = combined_vector_dec(Y_size + 1 : Y_size + Cb_size);
    zigzag_Cr_dec = combined_vector_dec(Y_size + Cb_size + 1 : Y_size + Cb_size+Cr_size);

    Y = izigzag(zigzag_Y_dec, m, n);
    Cb = izigzag(zigzag_Cb_dec, m, n);
    Cr = izigzag(zigzag_Cr_dec, m, n);
end