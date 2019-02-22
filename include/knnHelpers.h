#ifndef KNNHELPERS_H_4D8E793E
#define KNNHELPERS_H_4D8E793E


double doKnnClassification(const double *dist_mat_tril, const size_t n, const long long *true_class, const size_t k, const size_t n_perm, long long *knn_class, double *vote_matrix, double *per_perm);
ssize_t findUniqueItemsLL(const long long *list, const size_t n, long long *ulist);
size_t sub2ind_tril(size_t n_rows, size_t row, size_t col);
void Rpdist(const double* X, const size_t Rnx, const size_t Rp, double* distances);
void Rpdist_c(const double* X, const size_t n_rows, const size_t n_cols, double* distances);

#endif /* end of include guard: KNNHELPERS_H_4D8E793E */
