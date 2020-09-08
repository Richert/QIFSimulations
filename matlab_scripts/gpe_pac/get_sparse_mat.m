function C = get_sparse_mat(n1,n2,density,normalize)
%GET_SPARSE_MAT Summary of this function goes here
%   Detailed explanation goes here
C = randn(n1,n2);
vals_sorted = sort(reshape(C,n1*n2,1));
threshold = vals_sorted(int32(n1*n2*(1.0-density))+1);
C(C<threshold) = 0;
C(C>=threshold) = 1;
if normalize == 1
    row_sums = sum(C,2);
    for i=1:n1
        if row_sums(i) > 0
            C(i,:) = C(i,:) ./ row_sums(i);
        end
    end
end
end

