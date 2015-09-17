function p = multivariateGaussian(X, mu, Sigma2)
%MULTIVARIATEGAUSSIAN Computes the probability density function of the
%multivariate gaussian distribution.
%    p = MULTIVARIATEGAUSSIAN(X, mu, Sigma2) Computes the probability 
%    density function of the examples X under the multivariate gaussian 
%    distribution with parameters mu and Sigma2. If Sigma2 is a matrix, it is
%    treated as the covariance matrix. If Sigma2 is a vector, it is treated
%    as the \sigma^2 values of the variances in each dimension (a diagonal
%    covariance matrix)
%

k = length(mu);

if (size(Sigma2, 2) == 1) || (size(Sigma2, 1) == 1)
    Sigma2 = diag(Sigma2);
end
% 函数 bsxfun
% 【功能描述】两个数组间元素逐个计算.
% 【应用场合】当我们想对一个矩阵A的每一列或者每一行与同一个长度相等的向量a进行某些操作（比较大小，乘除等）时，我们只能用循环方法或者利用repmat函数将要操作的向量a复制成和A一样尺寸的矩阵，进而进行操作。从MATLAB R2007a开始，再遇到类似的问题时，我们有了简洁高效的方法，即利用bsxfun函数。
% 【函数描述】C=bsxfun(fun,A,B)：两个数组间元素逐个计算，fun是函数句柄或者m文件，也可以为如下内置函数 
% @plus 加 
% @minus 减 
% @times 数组乘 
% @rdivide 左除 
% @ldivide 右除 
X = bsxfun(@minus, X, mu(:)');
p = (2 * pi) ^ (- k / 2) * det(Sigma2) ^ (-0.5) * ...
    exp(-0.5 * sum(bsxfun(@times, X * pinv(Sigma2), X), 2));

end