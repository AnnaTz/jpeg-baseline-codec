function block = iBlockDCT(dctBlock)

% apply inverse dct to the given dequantized block
block = idct2(dctBlock);

end

