% do image fft
src = '../samples/1.jpg';
row = 3;
col = 2;
im = imread(src);
subplot(row, col, 1);
imshow(im);
gray = rgb2gray(im);
subplot(row, col, 2);
imshow(gray)
F = fft2(gray);
subplot(row, col, 3);
imshow(F, []);
FS = fftshift(F);
S = log(1 + abs(FS));
subplot(row, col, 4);
imshow(S, []);

fr = real(ifft2(ifftshift(FS)));
ret = im2uint8(mat2gray(fr));
subplot(row, col, 5);
imshow(ret);
