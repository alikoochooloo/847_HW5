format compact

% load data
load USPS.mat

% center the matrix
X = centering(A);

% calculate eigenvectors
V = eigenvectors(X);

% I used the following chcunck of code to produce the 2 subplots of the
% images reconstructed using pca
%{
V2 = V(:,1:10);
A2 = X*V2*V2';
A2 = reshape(A2(2,:), 16, 16);
subplot(1,5,1)
imshow(A2')
V2 = V(:,1:50);
A2 = X*V2*V2';
A2 = reshape(A2(2,:), 16, 16);
subplot(1,5,2)
imshow(A2')
V2 = V(:,1:100);
A2 = X*V2*V2';
A2 = reshape(A2(2,:), 16, 16);
subplot(1,5,3)
imshow(A2')
V2 = V(:,1:200);
A2 = X*V2*V2';
A2 = reshape(A2(2,:), 16, 16);
subplot(1,5,4)
imshow(A2')
A2 = reshape(A(2,:), 16, 16);
subplot(1,5,5)
imshow(A2')
%}

p = [10,50,100,200,256];
errors = [];

% do PCA baced on the number of eigenvectors used and record the errors
for i = 1:size(p,2)
    E = doNorm(V(:,1:p(i)),X);
    errors(i,:) = E;
end

% plot the resulting errors
plot(p,errors)

% function for centering the matrix columns
function X = centering(A)
    X = A;
    for j = 1:size(A,2)
        X(:,j) = A(:,j)-mean(A(:,j));
    end
end

% calculating and returning the eigenvectors
function V = eigenvectors(X)
    [U,S,V] = svd(X);
    
end

% doing Frobenius norm
function e = doNorm(V,X)
    do = X*V;
    undo = do * V';
    e = norm((X-undo), 'fro')^2;
end