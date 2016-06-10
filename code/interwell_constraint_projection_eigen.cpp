QList<Eigen::Vector3d> WellConstraintProjections::interwell_constraint_projection_eigen(QList<Eigen::Vector3d> coords, double d)
{
    /* If the two line segments already satisfy
     * the interwell distance constraint,
     * simply return the same coordinates
     */
    if( WellConstraintProjections::shortest_distance_eigen(coords) >=d){
        std::cout <<"Initial points satisfy constraints" << std::endl;
        return coords;
    }

    QList<Eigen::Vector3d> solution_coords;
    QList<Eigen::Vector3d> moved_coords;
    QList<Eigen::Vector3d> temp_coords;
    /* Iterate through moving points. First try moving 2 points, then 3 points
     * then 4 points. If problem can be solved moving k points, moving k+1 points
     * will be a worse solution. Return the best k point solution.
     */

    double cost = INFINITY;

    // Try moving 2 points
    std::cout <<"Initial points not feasible. Try moving 2 points" << std::endl;
    int two_point_index[4][2] = { {0, 2} ,
                                  {0, 3} ,
                                  {1, 2} ,
                                  {1, 3}  };

    for (int ii=0; ii<4; ii++){
        moved_coords = coords;
        temp_coords = WellConstraintProjections::well_length_projection_eigen(coords.at(two_point_index[ii][0]),coords.at(two_point_index[ii][1]),INFINITY,d,10e-5);
        moved_coords.replace(two_point_index[ii][0],temp_coords.at(0));
        moved_coords.replace(two_point_index[ii][1],temp_coords.at(1));
        if(WellConstraintProjections::shortest_distance_eigen(moved_coords) >=d && WellConstraintProjections::movement_cost_eig(coords,moved_coords) < cost){
            // If several moves of two points work, save the one with lovest movement cost
            cost = WellConstraintProjections::movement_cost_eig(coords,moved_coords);
            solution_coords = moved_coords;
        }
    }
    // If there were any succesful configurations, return the best one.
    if (cost < INFINITY){
        std::cout <<"Found 2-point solution" << std::endl;
        return solution_coords;
    }

    // ################## 3 POINT PART ############################
    // If no 2 point movements were succesful, try moving 3 points.
    std::cout <<"No 2 point solution. Try moving 3 points" << std::endl;
    int three_point_index[4][3] = { {2, 0, 1} ,
                                    {3, 0, 1} ,
                                    {0, 2, 3} ,
                                    {1, 2, 3}  };
    for (int ii=0; ii<4; ii++){
        // Reset moved coords to initial state
        moved_coords = coords;

        // Choose which 3 points to move. (order is important, second and third entry should belong to same line segment)
        QList<Eigen::Vector3d> input_cords_3p;
        for (int jj=0; jj<3; jj++){
            input_cords_3p.append(coords.at(three_point_index[ii][jj]));
        }
        Eigen::Matrix3d temp_A = WellConstraintProjections::build_A_3p_eigen(input_cords_3p);
        Eigen::Vector3d temp_b = WellConstraintProjections::build_b_3p_eigen(input_cords_3p,d);

        /* The kkt_eg_solutions solver handles some numerical issues
         * like A having some values close to machine epsilon and
         * eigenvalues being close to 0. Just assume that any solution
         * must be among the ones given in solution candidates. we check
         * all of them.
         */
        QList<Eigen::Vector3d> solution_candidates = WellConstraintProjections::kkt_eq_solutions_eigen(temp_A,temp_b);
        //std::cout << "there are " << solution_candidates.length() << " solution candidates" << std::endl;

        for(int sol_num = 0; sol_num < solution_candidates.length(); sol_num++){
            // Solution of three point problem
            temp_coords = WellConstraintProjections::move_points_3p_eigen(input_cords_3p,d,solution_candidates.at(sol_num));
            if(temp_coords.length()<1){temp_coords = input_cords_3p;}

            for (int jj=0; jj<3; jj++){
                moved_coords.replace(three_point_index[ii][jj],temp_coords.at(jj));
            }

            /*std::cout << "shortest distance unmoved 3p = " << shortest_distance_3p_eigen(input_cords_3p) << std::endl;
            std::cout << "shortest distance 3p = " << shortest_distance_3p_eigen(temp_coords) << std::endl;
            std::cout << "shortest distance 4p = " << shortest_distance_eigen(moved_coords) << std::endl;
            std::cout << "movement cost = " << WellConstraintProjections::movement_cost_eig(coords,moved_coords) << std::endl;*/

            if(WellConstraintProjections::shortest_distance_eigen(moved_coords) >=d-0.001 &&
                    WellConstraintProjections::movement_cost_eig(coords,moved_coords) < cost){
                // If several moves of two points work, save the one with lovest movement cost
                cost = WellConstraintProjections::movement_cost_eig(coords,moved_coords);
                solution_coords = moved_coords;
            }
        }

    }
    // If there were any succesful configurations, return the best one.
    if (cost < INFINITY){
        std::cout <<"Found 3-point solution" << std::endl;
        return solution_coords;}
    // ################## END 3 POINT PART ########################



    // ################## 4 POINT PART ############################
    std::cout <<"Found no 3-point solution. Try 4 points" << std::endl;

    // Get all candidates for vector s
    /* The kkt_eg_solutions solver handles some numerical issues
     * like A having some values close to machine epsilon and
     * eigenvalues being close to 0. Just assume that any solution
     * must be among the ones given in solution candidates. we check
     * all of them.
     */
    Eigen::Matrix3d temp_A = WellConstraintProjections::build_A_4p_eigen(coords);
    Eigen::Vector3d temp_b = WellConstraintProjections::build_b_4p_eigen(coords,d);
    QList<Eigen::Vector3d> solution_candidates = WellConstraintProjections::kkt_eq_solutions_eigen(temp_A,temp_b);

    // Go through candidates s and pick the best one
    for(int sol_num = 0; sol_num < solution_candidates.length(); sol_num++){

        moved_coords = WellConstraintProjections::move_points_4p_eigen(coords,d,solution_candidates.at(sol_num));
        if(WellConstraintProjections::shortest_distance_eigen(moved_coords) >=d-0.001 && WellConstraintProjections::movement_cost_eig(coords,moved_coords) < cost){
            // If several candidates for s work, save the one with lovest movement cost
            cost = WellConstraintProjections::movement_cost_eig(coords,moved_coords);
            solution_coords = moved_coords;
        }
    }

    if(solution_coords.length()>0){std::cout <<"Found 4-point solution" << std::endl;}
    else std::cout <<"Found no solution to problem" << std::endl;

    return solution_coords;
}