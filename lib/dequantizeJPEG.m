function dctBlock = dequantizeJPEG(qBlock,qTable,qScale)

% dequantize the given quantized block
% the final quantization table is the product of the given quantization
% table times the given quantization factor
dctBlock = qBlock .* (qScale * qTable);

end

