N = 60;
s = 10;

T = floor(N/s);
R = mod(N,s);

% create periodic basis vectors with period s
% these are the columns of the matrix X
I = eye(s);
X = repmat(I,T,1);
X = [X;I(1:R,:)];

% normalize X so columns have unit norm
X = X./sqrt(sum(X));

% look at Fourier transform
Y = fft(X)/sqrt(N);
