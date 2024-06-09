
This project was developed for the 'Multimedia Systems and Virtual Reality' course, part of the MEng in Electrical and Computer Engineering at the Aristotle University of Thessaloniki in 2019. The aim of the project is to develop a JPEG encoder/decoder following the ISO/IEC 10918-1:1994 standard, focusing on the baseline sequential DCT-based process. 

## Structure and Deliverables

### JPEG Library

The primary goal is to assemble a library of functions that encapsulate key components of the JPEG encoding and decoding template. Each function will include its inverse to facilitate reconstruction, integral to the decoding process.

#### Preprocessing

Initial preprocessing converts images from RGB to YCbCr format, supporting different subsampling formats (4:4:4, 4:2:2, and 4:2:0). The function signature is as follows:

```matlab
function [imageY, imageCb, imageCr] = convert2ycbcr(imageRGB, subsampling)
```

- **imageY, imageCb, imageCr**: Components of the YCbCr image.
- **imageRGB**: The original RGB image.
- **subsampling**: A 1x3 array defining the subsampling format (e.g., [4 2 0]).

The inverse function, responsible for converting YCbCr back to RGB and oversampling to match the original image's dimensions, is specified as:

```matlab
function imageRGB = convert2rgb(imageY, imageCr, imageCb, subsampling)
```

#### DCT Transformation

DCT transformations are applied at the block level:

```matlab
function dctBlock = blockDCT(block)
```

And its inverse:

```matlab
function block = iBlockDCT(dctBlock)
```

#### Quantization

Quantization of DCT blocks and the corresponding inverse function are defined as:

```matlab
function qBlock = quantizeJPEG(dctBlock, qTable, qScale)
function dctBlock = dequantizeJPEG(qBlock, qTable, qScale)
```

#### Zig-zag Scanning and RLE

Functions for run-length encoding (RLE) and its inverse are outlined for processing quantized DCT coefficients:

```matlab
function runSymbols = runLength(qBlock, DCpred)
function qBlock = irunLength(runSymbols, DCpred)
```

#### Huffman Coding

Huffman coding and decoding functions are provided for compressing run-length encoded symbols:

```matlab
function huffStream = huffEnc(runSymbols)
function runSymbols = huffDec(huffStream)
```

### Demo 1

`demo1.m` demonstrates the process of converting an RGB image to YCbCr and back, along with DCT coefficient processing. Different subsampling formats are used for two provided images.

## JPEG Integration

Integrating the developed functions into a comprehensive encoder/decoder function is the second major deliverable, aiming at both quantitative and qualitative analysis of the achieved compression.

### JPEG Encoder/Decoder

The function signatures for encoding and decoding are:

```matlab
function JPEGenc = JPEGencode(img, subsampling, qScale)
function imgRec = JPEGdecode(JPEGenc)
```

### Demo 2

`demo2.m` will calculate and compare entropies in different domains for two provided images, using the same parameters as in Demo 1.

## JPEG Syntax

The final deliverable focuses on constructing and decoding a bitstream of the encoded image, ensuring compatibility with standard image processing software.

### JPEG Syntax Encoder/Decoder

The functions for encoding to and decoding from a JPEG bitstream are defined as:

```matlab
function JPEGencStream = JPEGencodeStream(img, subsampling, qScale)
function imgCmp = JPEGdecodeStream(JPEGencStream)
```
