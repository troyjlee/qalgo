N = 21^2;
% s should be the order of x for the x you choose
s = 6;

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

% sampling
% sample a random coset state
j = randi(s);

% sample entry b with prob |Y(b,j)|^2
probs = abs(Y(:,j)).^2;
thresh = rand(1);
% create vector of partial sums of probs
sum_probs = tril(ones(N))*probs;
% find the first b with sum_probs(b) > thresh
[arr,indices] = sort([thresh;sum_probs]);
disp('The measurement gives index');
b = find(indices == 1)

