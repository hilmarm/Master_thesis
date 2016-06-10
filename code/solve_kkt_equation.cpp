QList<Eigen::Vector3d> kkt_eq_solutions_eigen(Eigen::Matrix3d A, Eigen::Vector3d b)
{
    
    QList<Eigen::Vector3d> candidate_solutions;

    /* First assume that A-\mu I has an inverse.
     * We can find the inverse of it and solve
     * a sixth degree equation for \mu.
     */

    // Remove values that are close to zero. Mostly caused by initialization.
    A = WellIndexCalculator::WellConstraintProjections::rm_entries_eps_matrix(A,10e-12);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> A_es(A);

    // Need to remove eigenvalues which are close to zero
    Eigen::Vector3d eigenvalues = WellConstraintProjections::rm_entries_eps(A_es.eigenvalues(),10e-12);

    // Calculate the coefficients of the sixth degree polynomial
    Eigen::VectorXd coeffs = WellConstraintProjections::coeff_vector_eigen
                                    (eigenvalues, A_es.eigenvectors().inverse(), b);

    /* There is an issue where coefficients should be zero but are not
     * but because of numerical issues these need to be handled manually.
     * Simply set all whose fabs(x)<10-e12 to zero.
     */
    coeffs = WellConstraintProjections::rm_entries_eps_coeffs(coeffs,10e-12);


    // Compute roots of polynomial with RPOLY library
    Eigen::VectorXd realroots(6);
    Eigen::VectorXd comproots(6);
    rpoly_plus_plus::FindPolynomialRootsJenkinsTraub(coeffs,&realroots,&comproots);

    // Loop through all roots of polynomial
    for (int ii=0;ii<6;ii++){

        // Root may not be complex or an eigenvalue of A
        if(comproots[ii]==0 && eigenvalues[0]!=realroots[ii] &&
                eigenvalues[1]!=realroots[ii] && eigenvalues[2]!=realroots[ii]){

            // We have found a valid root. Get vector s.
            double cur_root = realroots[ii];
            Eigen::Vector3d cur_root_vec;
            cur_root_vec << cur_root, cur_root, cur_root;
            Eigen::Matrix3d invmatr = (eigenvalues-cur_root_vec).asDiagonal();

            // Get vector candidate vector s
            Eigen::Vector3d s = A_es.eigenvectors()*invmatr.inverse()*A_es.eigenvectors().inverse()*b;
            candidate_solutions.append(s);
        }
    }

    /* Now for the second part assume that A-\mu I is not
     * invertible, i.e. \mu is an eigenvalue of A. Then
     * we either have an infinite amount of solutions of
     * (A-\mu I)s = b. Require s have length 1 to find
     * at most two solutions as long as all points are
     * not on the same line.
     */

    // Loop through all 3 eigenvalues of A
    for(int i=0; i<3; i++){

        QList<Eigen::Vector3d> eigenvalue_solutions;

        // Create linear system (A-\my I)s = b
        Eigen::Matrix3d A_eig = A- eigenvalues[i]*Eigen::Matrix3d::Identity();
        Eigen::Vector3d b_eig = b;

        // Check for existence of solutions
        if(WellConstraintProjections::solution_existence(A_eig,b_eig)){

            // Solves non-invertible case and returns, if any, the feasible vectors s
            eigenvalue_solutions = WellConstraintProjections::non_inv_solution(A_eig, b_eig);
        }

        // If any solutions exist, add them to solution_vectors
        for(int jj = 0; jj < eigenvalue_solutions.length(); jj++){
            candidate_solutions.append(eigenvalue_solutions.at(jj));
        }
    }

    return candidate_solutions;
}