format compact

% creating random samples
x = [rand(1,50)*2+3.1,rand(1,50)*2+4.4];
y = [rand(1,50)*2+3.1,rand(1,50)*2+4.4];
A = [x',y'];

%calling the clustering function
[C,M] = clustering(A,2);

%distributing each sample according to their cluster
[D1,D2] = redistribution(A,M);

%ploting the samples and the centers
plot(D1(:,1),D1(:,2),'.','MarkerSize', 16)
hold on
plot(D2(:,1),D2(:,2),'.','MarkerSize', 16)
hold on
plot(C(:,1),C(:,2),'.','MarkerSize', 26)
hold off

% redistribution function creates 2 arrays for the two clusters
function [D1,D2] = redistribution(A,M)
    D1 = [];
    D2 = [];
    for i = 1:size(A,1)
        if M(i,1)
            D1(size(D1,1)+1,:) = A(i,:);
        else
            D2(size(D2,1)+1,:) = A(i,:);
        end
    end
end

% main clustering function that calls for M and C to be updated in 6
% itterations and diplay the error and spectral relaxation error of the matrix
function [c,m] = clustering(A,k)
    cx = rand(k,1)*3.5+3.1;
    cy = rand(k,1)*3.5+3.1;
    c = [cx,cy];
    for i = 1:6
        m = calcM(A,c,k);
        c = calcC(A,m,k);
        disp(standarderror(A,m,c,k));
    end
    disp(spectralrelaxation(A));
end

% spectral relaxation calculation error function
function error = spectralrelaxation(A)
    [U,S,V] = svd(A);
    j = min(size(A));
    sum = 0;
    for i = 1:j
        sum = sum + S(i,i)^2;
    end
    tr = trace(A'*A);
    error = tr - sum;
end

% the standart error calculation function
function sums = standarderror(A,M,C,K)
    sums = 0;
    for j = 1:K
        for i = 1:size(M,1)
            if M(i,j)
                sub = A(i,:)-C(j,:);
                sums = sums + (norm(sub))^2;
            end
            
        end
    end
end

% calculating the M matrix
function M =  calcM(A,C,k)
    M = zeros(size(A,1),k);
    for i = 1:size(A,1)
        ind = closerC(A(i,:),C,k);
        M(i,ind) = 1;
        
    end 
end

% finding which center is closest to the sample
function cluster = closerC(a,C,k)
    cluster = 0;
    measure = 1/0;
    for j = 1:k
        sub = a-C(j,:);
        
        localmeasure = norm(sub)^2;
        if measure > localmeasure
            cluster = j;
            measure = localmeasure;
        end
    end
end

% calculating the cluster centers
function C = calcC(A,M,K)
    C = zeros(K,size(A,2));
    for j = 1:K
        sums = zeros(1,size(A,2));
        count = 0;
        for i = 1:size(A,1)
            if M(i,j)
                sums = sums+A(i,:);
                count = count +1; 
            end
        end
        
        C(j,:) = sums/count;
    end
end
