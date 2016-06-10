QList<Eigen::Vector3d> WellConstraintProjections::well_length_projection_eigen(
    Eigen::Vector3d heel, Eigen::Vector3d toe, double max, double min, double epsilon)
{
    QList<Eigen::Vector3d> projected_points;
    Eigen::Vector3d moved_heel;
    Eigen::Vector3d moved_toe;

    // Need the vector going from heel to toe to project points
    Eigen::Vector3d heel_to_toe_vec = toe-heel;
    // Distance between heel and toe.
    double d = heel_to_toe_vec.norm();

    // heel and toe same point.
    // All directions equally good.
    if(d==0){
        Eigen::Vector3d unit_vector;
        unit_vector << 1,0,0;
        moved_heel = heel + (min/2)*unit_vector;
        moved_toe  = heel - (min/2)*unit_vector;
        projected_points.append(moved_heel);
        projected_points.append(moved_toe);
        return projected_points;
    }

    // Normalize vector to get correct distance
    heel_to_toe_vec.normalize();

    // Trivial case
    if(d<=max && d>=min){
        projected_points.append(heel);
        projected_points.append(toe);
    }

    // Distance too long
    else if(d>max){
        double move_distance = 0.5*(d-max+(epsilon/2));
        moved_heel = heel + move_distance*heel_to_toe_vec;
        moved_toe =  toe  - move_distance*heel_to_toe_vec;
        projected_points.append(moved_heel);
        projected_points.append(moved_toe);
    }

    // Distance too short
    else{
        double move_distance = 0.5*(d-min-(epsilon/2));
        moved_heel = heel + move_distance*heel_to_toe_vec;
        moved_toe =  toe  - move_distance*heel_to_toe_vec;
        projected_points.append(moved_heel);
        projected_points.append(moved_toe);
    }

    return projected_points;
}