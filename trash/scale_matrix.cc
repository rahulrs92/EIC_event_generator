#include<detector_properties.hh>

template<size_t M, size_t N>
int detector_properties::scale_matrix(const matrix<M, N>& Y, matrix<M, N>& X, matrix<M, 1>& alpha, matrix<N, 1>& beta)
{
    constexpr size_t K = M + N - 1; // 7 parameters (b_1 fixed to 0 for gauge freedom)

    // Target constraints
    matrix<M, 1> r = N;
    matrix<N, 1> c = M;

    // Log-space parameters vector p = [a_0, ..., a_5, b_0]^T
    matrix<K, 1> p(0.0);

    double tol = 1e-10;
    int max_iter = 50;
    int iter = 0;

    for (; iter < max_iter; ++iter)
    {
        // 1. Calculate current scaled yields: X_ik = Y_ik * exp(a_i + b_k)
        for (size_t i = 0; i < M; ++i)
        {
            X(i, 0) = Y(i, 0) * std::exp(p[i] + p[M]); // b_0 = p[M]
            X(i, 1) = Y(i, 1) * std::exp(p[i]);        // b_1 = 0
        }

        // 2. Compute Row (R) and Column (C) sums
        matrix<M, 1> R(0.0);
        matrix<N, 1> C(0.0);
        for (size_t i = 0; i < M; ++i)
        {
            R[i] = X(i, 0) + X(i, 1);
            C[0] += X(i, 0);
            C[1] += X(i, 1);
        }

        // 3. Construct residual gradient vector g
        matrix<K, 1> g(0.0);
        double max_err = 0.0;
        for (size_t i = 0; i < M; ++i)
        {
            g[i] = R[i] - r[i];
            max_err = std::max(max_err, std::abs(g[i]));
        }
        g[M] = C[0] - c[0];
        max_err = std::max(max_err, std::abs(g[M]));

        if (max_err < tol) break;

        // 4. Build Hessian matrix H
        matrix<K, K> H(0.0);

        // Diag(R) and off-diagonal coupling with column 0
        for (size_t i = 0; i < M; ++i)
        {
            H(i, i) = R[i];
            H(i, M) = X(i, 0);
            H(M, i) = X(i, 0);
        }
        // Diag(C_sub)
        H(M, M) = C[0];

        // 5. Compute Newton step using matrix division operator (solves H * dp = -g)
        matrix<K, 1> dp = (-g) / H;

        // 6. Update parameters
        p += dp;
    }

    // Extract multipliers
    for (size_t i = 0; i < M; ++i) alpha[i] = std::exp(p[i]);
    beta[0] = std::exp(p[M]);
    beta[1] = 1.0;
    return iter;
}
