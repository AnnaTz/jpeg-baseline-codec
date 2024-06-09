function qBlock = quantizeJPEG(dctBlock,qTable,qScale)

% quantize the given dct block
% the final quantization table is the product of the given quantization
% table times the given quantization factor
qBlock = dctBlock ./ (qScale * qTable);

end

