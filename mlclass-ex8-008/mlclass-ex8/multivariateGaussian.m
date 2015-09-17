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
% ���� bsxfun
% ���������������������Ԫ���������.
% ��Ӧ�ó��ϡ����������һ������A��ÿһ�л���ÿһ����ͬһ��������ȵ�����a����ĳЩ�������Ƚϴ�С���˳��ȣ�ʱ������ֻ����ѭ��������������repmat������Ҫ����������a���Ƴɺ�Aһ���ߴ�ľ��󣬽������в�������MATLAB R2007a��ʼ�����������Ƶ�����ʱ���������˼���Ч�ķ�����������bsxfun������
% ������������C=bsxfun(fun,A,B)�����������Ԫ��������㣬fun�Ǻ����������m�ļ���Ҳ����Ϊ�������ú��� 
% @plus �� 
% @minus �� 
% @times ����� 
% @rdivide ��� 
% @ldivide �ҳ� 
X = bsxfun(@minus, X, mu(:)');
p = (2 * pi) ^ (- k / 2) * det(Sigma2) ^ (-0.5) * ...
    exp(-0.5 * sum(bsxfun(@times, X * pinv(Sigma2), X), 2));

end