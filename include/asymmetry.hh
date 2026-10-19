#include<matrix3.hh>


namespace Polarimeter
{
    /**
     * @brief Struct to hold extracted flux, individual asymmetries, and full covariance matrix.
     * 
     * @tparam M Number of azimuthal detectors in the polarimeter array.
     */
    template<size_t M>
    struct PolarimeterMetrics
    {
        double L;               // Relative beam flux ratio: L = lambda_1 / lambda_2
        double L_e;             // 1-sigma uncertainty on relative flux ratio (sigma_L)
        matrix<M, 1> A;         // Column vector of exact detector asymmetries A_i (M x 1)
        matrix<M, M> cov_A;     // Full M x M covariance matrix cov(A_i, A_j) including flux correlations
        matrix<2, 2> cov_Ab;    // 2 x 2 covariance matrix cov(A_bs, A_bc) of physical beam parameters
        int status = 0;         // Status code: 0 = success, -1 = non-positive counts/errors
    };

    // size_t constexpr M = 12; // Number of detectors in the polarimeter array

    /**
     * @brief Adaptively computes relative flux L, individual asymmetries A_i, covariance matrix cov_A,
     *        and physical beam polarization parameters Ab = [A_bs, A_bc]^T.
     * 
     * METHODOLOGY:
     * 1. Relative Flux Estimation (L):
     *    - Even M: Closed-form geometric mean ratio across antipodal detector pairs.
     *    - Odd M >= 3: Weighted log-ratio linear regression fitting ln(y_1i / y_2i) = ln L + 2 A_bs sin(phi_i) + 2 A_bc cos(phi_i).
     * 2. Exact Detector Asymmetries A_i:
     *    - Calculated directly as A_i = (y_1i - L*y_2i) / (y_1i + L*y_2i).
     * 3. Exact First-Order Error Propagation (cov_A):
     *    - Formulates total Jacobian J = J_direct + beta * J_L^T for stacked yield vector Y = [y_1; y_2].
     *    - Computes full M x M asymmetry covariance matrix C = J * Sig_Y * J^T.
     * 4. Generalized Least Squares (GLS) Fit (Ab):
     *    - Inverts full covariance matrix C to form precision matrix W = C^-1.
     *    - Solves physical model A_i = A_bs sin(phi_i) + A_bc cos(phi_i) taking all inter-detector correlations into account.
     * 
     * @tparam M Number of detectors (M >= 2)
     * @param Y   Input M x 2 yield counts: Col 0 = Beam UP (y_1i), Col 1 = Beam DOWN (y_2i)
     * @param Y_e Input M x 2 measurement uncertainties (s_1i, s_2i)
     * @param R   Input M x 2 Geometric Design matrix: Col 0 = Ay*sin(phi_i), Col 1 = Ay*cos(phi_i)
     * @param Ab  Output 2 x 1 vector of extracted physical beam asymmetries [A_bs, A_bc]^T
     * @param Ab_e Output 2 x 1 vector of parameter uncertainties [sigma(A_bs), sigma(A_bc)]^T
     * @return PolarimeterMetrics<M> Struct containing L, sigma_L, A_i, cov(A_i), cov(Ab), and status.
     */
    template<size_t M>
    PolarimeterMetrics<M> scale_matrix
    // PolarimeterMetrics<M> scale_matrix
    (
        const matrix<M, 2>& Y, 
        const matrix<M, 2>& Y_e, 
        const matrix<M, 2>& R, 
        matrix<2, 1>& Ab,
        matrix<2, 1>& Ab_e
    ) 
    {
        static_assert(M >= 2, "Error: At least 2 detectors are required.");

        PolarimeterMetrics<M> res;

        // Validate inputs: Raw yields and uncertainties must be strictly positive
        for (size_t i = 0; i < M; ++i)
        {
            if (Y(i,0) <= 0.0 || Y(i,1) <= 0.0 || Y_e(i,0) <= 0.0 || Y_e(i,1) <= 0.0)
            {
                res.status = -1;
                return res;
            }
        }

        const double pi = M_PI;

        // Aliases for clear output variable references
        auto &L   = res.L;
        auto &L_e = res.L_e;
        auto &C   = res.cov_A;
        auto &A   = res.A;
        auto &CAb = res.cov_Ab;
        
        // Initialize outputs
        L    = 0.0;
        L_e  = 0.0;
        Ab   = 0.0;
        Ab_e = 0.0;

        // =====================================================================
        // STEP 1: Compute Relative Beam Flux L and its Variance sigma_L^2
        // =====================================================================

        matrix<M,3> Mlog;     // Log-ratio regression design matrix: [1, 2*sin(phi_i), 2*cos(phi_i)]
        matrix<M,M> Wlog;     // Diagonal weight matrix for log ratios: w_i = 1 / ( (s_1i/y_1i)^2 + (s_2i/y_2i)^2 )
        matrix<3,3> Clog;     // Log-ratio parameter covariance matrix: (Mlog^T * Wlog * Mlog)^-1
        matrix<M,1> dlog;     // Log-ratio measurement vector: d_i = ln(y_1i / y_2i)
        matrix<M,1> U0log;    // Projection weight vector extracting L sensitivity: U_0 = Wlog * Mlog * Clog * [1,0,0]^T
        matrix<3,1> plog;     // Log-ratio fitted parameters: [ln L, A_bs, A_bc]^T

        if (M % 2 == 0)
        {
            // -----------------------------------------------------------------
            // EVEN M: Closed-form initial guess via geometric mean
            // -----------------------------------------------------------------
            for (size_t i = 0; i < M; ++i)
            {
                L += log(Y(i,0) / Y(i,1));
                
                // Accumulate relative variance terms: sum_i [ (s_1i/y_1i)^2 + (s_2i/y_2i)^2 ]
                for (size_t j = 0; j < 2; ++j) L_e += (Y_e(i,j) * Y_e(i,j)) / (Y(i,j) * Y(i,j));
            }
            L = exp(L / M);                                  // Final relative flux L

            // -----------------------------------------------------------------
            // Newton-Raphson Refinement: solves sum_i A_i(L) = 0
            // -----------------------------------------------------------------
            // for (size_t iter = 0; iter < 3; ++iter)
            // {
            //     double a = 0.0, b = 0.0;
            //     for (size_t i = 0; i < M; ++i)
            //     {
            //         double Di = Y(i,0) + L * Y(i,1);
            //         a += (Y(i,0) - L * Y(i,1)) / Di;
            //         b += 2.0 * Y(i,0) * Y(i,1) / (Di * Di);
            //     }
            //     L += a / b;
            // }

            // -----------------------------------------------------------------
            // Final uncertainty L_e using refined L
            // -----------------------------------------------------------------
            L_e = sqrt((L * L / (M * M)) * L_e);                // sigma_L = (L / M) * sqrt( sum_i relative_variances )
        }
        else
        {
            // -----------------------------------------------------------------
            // ODD M >= 3: Log-ratio linear regression to extract ln(L)
            // System equation: ln(y_1i / y_2i) = ln L + 2*A_bs*sin(phi_i) + 2*A_bc*cos(phi_i)
            // -----------------------------------------------------------------
            for (size_t i = 0; i < M; ++i)
            {
                double phi_i = 2.0 * pi * static_cast<double>(i) / static_cast<double>(M);
                
                Mlog[i] = {1.0, 2.0 * R(i,0), 2.0 * R(i,1)};  // Design matrix row: [1, 2*sin(phi_i), 2*cos(phi_i)]
                
                // Weight w_i is inverse variance of ln(y_1i / y_2i) via log-variance transform
                Wlog[i][i] = 1.0 / (pow(Y_e(i,0)/Y(i,0), 2) + pow(Y_e(i,1)/Y(i,1), 2));
                
                dlog[i] = log(Y(i,0) / Y(i,1));
            }

            // Fit parameters: plog = (Mlog^T * Wlog * Mlog)^-1 * Mlog^T * Wlog * dlog
            Clog = (Mlog.tp() * Wlog * Mlog) ^ -1;
            plog = Clog * (Mlog.tp() * Wlog * dlog);

            L   = exp(plog[0]);                 // plog[0] = ln L  =>  L = exp(plog[0])
            L_e = L * sqrt(Clog(0,0));          // Error propagation: sigma_L = L * sigma_(ln L)

            // Matrix projection vector U0log = Wlog * Mlog * Clog * [1, 0, 0]^T
            // Encapsulates sensitivity of ln L to each individual detector log-ratio dlog[i]
            matrix<3,1> tr = {1.0, 0.0, 0.0};
            U0log = Wlog * Mlog * Clog * tr;
        }
        
        // =====================================================================
        // STEP 2: Compute Exact Physical Asymmetries A_i
        // Formula: A_i = (y_1i - L * y_2i) / (y_1i + L * y_2i)
        // =====================================================================
        
        // Di is an M x 1 column vector of denominators: D_i = y_1i + L * y_2i
        auto Di = Y(0) + L * Y(1);

        for (size_t i = 0; i < M; ++i) A[i] = (Y(i,0) - L * Y(i,1)) / Di[i];

        // =====================================================================
        // STEP 3: Full Asymmetry Covariance Matrix C = cov(A_i, A_j)
        // First-Order Error Propagation: C = J * Sig_Y * J^T
        // =====================================================================
        
        matrix<M, 2*M> J;         // Full Jacobian matrix of A_i w.r.t stacked yield vector Y = [y_1; y_2]
        matrix<2*M, 2*M> Sig_Y;   // Input yield covariance matrix (2M x 2M diagonal)
        matrix<M, 1> beta;        // Vector of flux sensitivity partial derivatives: beta_i = dA_i / dL
        matrix<M, 1> h;           // Vector for cross-correlation coupling terms
        matrix<2*M, 1> Jl;        // Gradient vector of beam flux w.r.t yields: Jl = [dL/dy_1; dL/dy_2]
        
        for (size_t i = 0; i < M; ++i)
        {
            // 1. Diagonal yield measurement variances
            Sig_Y(i, i)         = Y_e(i,0) * Y_e(i,0);    // var(y_1i) = s_1i^2
            Sig_Y(i + M, i + M) = Y_e(i,1) * Y_e(i,1);    // var(y_2i) = s_2i^2

            // 2. Direct yield derivatives (holding flux L constant)
            J(i, i)     =  2.0 * L * Y(i,1) / pow(Di[i], 2);    // alpha_i = (dA_i / dy_1i)_direct
            J(i, i + M) = -2.0 * L * Y(i,0) / pow(Di[i], 2);    // gamma_i = (dA_i / dy_2i)_direct

            // 3. Direct flux derivative
            beta[i] = -2.0 * Y(i,0) * Y(i,1) / pow(Di[i], 2);   // beta_i = dA_i / dL

            // 4. Flux gradient terms dL / dy_1i and dL / dy_2i
            if (M % 2 == 0)
            {
                // Even M: Derived from geometric mean flux estimator
                Jl[i]     =  L / (static_cast<double>(M) * Y(i,0));    // dL / dy_1i
                Jl[i + M] = -L / (static_cast<double>(M) * Y(i,1));    // dL / dy_2i
            }
            else
            {
                // Odd M: Derived from log-ratio linear regression weights U0log
                Jl[i]     =  L * U0log[i] / Y(i,0);                    // dL / dy_1i
                Jl[i + M] = -L * U0log[i] / Y(i,1);                    // dL / dy_2i
            }
        }

        // 5. Total Jacobian assembly: J_total = J_direct + beta * Jl^T
        // Adds implicit flux variance and inter-detector correlation contributions
        J += beta * Jl.tp(); 

        // 6. Asymmetry Covariance Matrix Sandwich Product
        C = J * Sig_Y * J.tp(); 

        // =====================================================================
        // STEP 5: Final Generalized Least Squares (GLS) Fit for Beam Parameters
        // Fits physical model: A_i = A_bs * sin(phi_i) + A_bc * cos(phi_i)
        // =====================================================================
        
        // Construct physical design matrix Mphys (M x 2)
        matrix<M, 2> Mphys;
        for (size_t i = 0; i < M; ++i)  Mphys[i] = {R(i,0), R(i,1)};  // Col 0 = sin(phi_i), Col 1 = cos(phi_i)
        
        // 1. Precision Matrix W = C^-1 (Inverts full M x M covariance matrix including off-diagonals)
        // auto W = C ^ -1;
        matrix<M, M> W;
        double sf = C.norm() / M; // Scale factor for numerical stability
        // cout << "C = " << C << endl;
        
        // Covariance matrix is guaranteed to have at least one singularity due to the flux normalization constraint.
        // Full pseudo inversion using SVD to handle singular covariance matrices
        W = (C / sf).pseudo_inverse(1e-10); // Use pseudo-inverse with threshold for numerical stability
        W /= sf; // Rescale back to original units

        // Diagonal regularization fallback
        // for (size_t i = 0; i < M; ++i) W(i,i) = 1.0 / (C(i,i) + 1e-10);
        // cout << "Precision matrix W = C^-1 diagonal:\n" << W.diag() << endl;

        // Unweighted fallback (identity matrix)
        // W.diag() = 1.0 / sf; // Use scaled identity matrix as a fallback
        // cout << "Precision matrix W = C^-1 diagonal:\n" << W.diag() << endl;
        
        // 2. Physics Covariance Matrix (2 x 2)
        CAb = (Mphys.tp() * W * Mphys) ^ -1;
        
        // 3. Aitken's GLS Estimator (2 x 1)
        Ab = CAb * (Mphys.tp() * W * A);
        
        // 4. Store final physical parameter uncertainties
        Ab_e[0] = sqrt(CAb(0, 0));     // sigma(A_bs)
        Ab_e[1] = sqrt(CAb(1, 1));     // sigma(A_bc)

        res.status = 0;
        return res;
    }

    // template <size_t M>
    // matrix<M,M> C2W_singular(const matrix<M,M>& C)
    // {
    //     matrix<M,M> W;
    //     double sf = C.trace() / M; // Scale factor for numerical stability
    //     // Full pseudo inversion fallback using SVD for singular covariance matrices
    //     cout << "Warning: Singular covariance matrix detected. Using Singular Value Decomposition for precision matrix." << endl;
    //     matrix<M, M> U, iS, V;
    //     valarray<double> Sv;
    //     (C / sf).SVD(U, Sv, V);
    //     cout << "Singular values of C / sf:\n" << Sv << endl;
    //     for (size_t i = 0; i < M; ++i) iS(i,i) = abs(Sv[i]) > 1e-10 ? 1.0 / Sv[i] : 0.0;
    //     W = V * iS * U.tp();
    //     // valarray<double> diagW = W.diag();
    //     // diagW /= diagW.sum() / M; // Normalize diagonal for comparison
    //     // matrix<M,1> diagWm; for (size_t i = 0; i < M; ++i) diagWm[i] = diagW[i];
    //     // cout << "Efficiency estimate from weights:\n" << diagWm << endl;
    //     W /= sf;
    //     cout << "Precision matrix W = C^-1 diagonal:\n" << W.diag() << endl;
    //     return W;
    // }


    // Add at the bottom of asymmetry.cc:
    // template PolarimeterMetrics<N_R> detector_properties::scale_matrix<N_R>(
    //     const matrix<N_R, 2>&, 
    //     const matrix<N_R, 2>&, 
    //     const matrix<N_R, 2>&, 
    //     matrix<2, 1>&, 
    //     matrix<2, 1>&
    // );

} // namespace Polarimeter