% parameters: you can change these to play around
% We work over Z_M
% and compute log_g(x)
M = 23;
N = sum(gcd(1:M-1,M)==1);
g = 3;
x = 18;

% First we make the matrix for the  Fourier transform over Z_N
clear i;
omega = exp(2*pi*i/N);
elts = [0:N-1];

c = omega.^-elts;
F = fliplr(vander(c))/sqrt(N);

% check that it works
f = randn(N,1);
max(abs(F*f - fft(f)/sqrt(N)))

F2 = kron(F,F);

% discrete log portion

% compute discrete log in the dumb way
logx = -1;
for j=0:M-1
  if mod(g^j,M) == x
    logx = j;
    break
  endif
endfor
if logx == -1
  disp('no power of g equals x mod M')
  return
else
  disp('the answer should be')
  logx
endif

% build the subgroup
% H = [elts',mod(-logx*elts,N)'];

% Build a matrix X whose columns have the second 
% element of the pair of the cosets of H

X = (-logx*elts)' + elts;
X = mod(X,N);

% make the characteristic vectors of each column of X
% these will be N^2 dimensional vectors
% map group element (a,b) to a*N + b + 1
R = N*elts'*ones(1,N);
row_indices = 1+(R+X)(:);
col_indices = (ones(N,1)*[1:N])(:);
Y = sparse(row_indices,col_indices, (1/sqrt(N))*ones(N^2,1),N^2,N);

% Look at Fourier transform over Z_N^2 of the coset states
Yhat = F2*Y;

% sample a random coset state
j = randi(N);
probs = abs(Yhat(:,j)).^2;
thresh = rand(1);
% create vector of partial sums of probs
sum_probs = tril(ones(N^2))*probs;
% find the first b with sum_probs(b) > thresh
[arr,indices] = sort([thresh;sum_probs]);
index = find(indices == 1);
%Nits = dec2base(index-1,N,2) - '0';
Nits = zeros(2,1);
Nits(2) = mod(index-1,N);
Nits(1) = (index-1-Nits(2))/N;
[val,binv] = gcd(Nits(2),N);
if val == 1
  disp('success');
  binv = mod(binv,N);
  answer = mod(binv*Nits(1),N)
else
  disp('unlucky measurement');
endif
