#ifndef WELL_CONSTRAINT_PROJECTIONS_H
#define WELL_CONSTRAINT_PROJECTIONS_H


#include <RpolyPlusPlus/polynomial.h>
#include <RpolyPlusPlus/find_polynomial_roots_jenkins_traub.h>
#include "Model/reservoir/grid/xyzcoordinate.h"
#include "WellIndexCalculator/geometry_functions/geometryfunctions.h"
#include "WellIndexCalculator/geometry_functions/geometryfunctions_exceptions.h"
#include <QList>
#include <armadillo>
#include <Eigen/Dense>
#include <QVector3D>
#include <QList>

namespace WellIndexCalculator {
/*!
 * \brief WellConstraintProjections is a collection for solving projection of well constraints
 *
 * We want to collect all mathematical functions for projection wells back into feasible space
 * Most of these functions will be used for the interwell distance problem as it concerns itself
 * with finding and solving a sextic equation
 */
class WellConstraintProjections
{

public:
    // simple summation and product functions
    static double sum_i(double a,double b,double c);
    static double sum_ij(double a,double b,double c);
    static double prod_i(double a,double b,double c);
    static void real_roots(Eigen::VectorXd coeffs, Eigen::VectorXd *real_roots);

    // Computes the candidate vectors that solve interwell distance problem
    static QList<QVector3D*> compute_interwell_vector_4d(QList<Model::Reservoir::Grid::XYZCoordinate*> coords, double d);

    static QList<QVector3D *> interwell_projection_4p(QList<QVector3D*> coords, double d);

    // Computes the candidate vectors that solve interwell distance problem
    static QList<QVector3D*> compute_interwell_vector_3d(QList<Model::Reservoir::Grid::XYZCoordinate*> coords, double d);
    static QList<QVector3D *> interwell_projection_3p(QList<QVector3D*> coords, double d);
    static QList<QVector3D *> interwell_projection_3p_eigen(QList<QVector3D*> coords, double d);


    static QList<QVector3D*> well_distance_projection(QVector3D* heel, QVector3D* toe, double max, double min, double epsilon);

    static QList<QVector3D*> move_points_4p(QList<QVector3D*> coords, double d, QVector3D *s);
    static QList<QVector3D*> move_points_3p(QList<QVector3D*> coords, double d, QVector3D *s);

    static QList<Model::Reservoir::Grid::XYZCoordinate*> move_points_4p_XYZ(QList<Model::Reservoir::Grid::XYZCoordinate*> coords, double d, QVector3D *s);

    static QList<QVector3D *> kkt_eq_solutions(arma::mat A, arma::vec b);

    static double shortest_distance_3p(QList<QVector3D *> coords);

    /*!
     * \brief Master function for the interwell constraint projection. Projects two wells so that they are at least a distance d appart.
     * \param Two line segments defined by their endpoints. Two first entries belong to first line segment, two last entries belong to second line segment.
     * \return Cheapest(L2-norm) projection of line segments s.t. they are at least a distance d appart.
    */
    static QList<QVector3D*> interwell_constraint_projection(QList<QVector3D*> coords, double d);

    /*!
     * \brief Master function for the interwell constraint projection. Projects two wells so that they are at least a distance d appart.
     * \param Two line segments defined by their endpoints. Two first entries belong to first line segment, two last entries belong to second line segment.
     * \return Cheapest(L2-norm) projection of line segments s.t. they are at least a distance d appart.
    */
    static QList<QVector3D*> interwell_constraint_projection_2nd(QList<QVector3D*> coords, double d);

    /*!
     * \brief Calculates the shortest distance between two line segments. Line segment endpoints given by 4 coordinates
     * \param coordinates containing endpoints of both line segments
     * \return shortest distance between line segments
    */
    static double shortest_distance(QList<QVector3D*> coords);

    static QList<Eigen::Vector3d> non_inv_solution(Eigen::Matrix3d A, Eigen::Vector3d b);

    static bool solution_existence(Eigen::Matrix3d A, Eigen::Vector3d b);

    // Below here used to be private

    /*!
     * \brief returns the k'th element of the matrix-vector product Qt*b
     * \param matrix
     * \param vector
     * \param k'th element
     * \return k'th element of matrix-vector product
    */
    static double prod_Qtb(arma::mat Qt, arma::vec b, int k);

    /*!
     * \brief Computes eigenvectors of a symmetric matrix and collects them in a matrix which is transposed before returning it
     * \param matrix
     * \return transposed matrix of eigenvectors
    */
    static arma::mat Q_trans(arma::mat A);
    // Build matrices for 2point and 3point case!

    static arma::mat build_A(QList<QVector3D*> coords);
    static arma::vec build_b(QList<QVector3D*> coords, double d);

    static arma::mat build_A_4p_qvec(QList<QVector3D*> coords);
    static arma::vec build_b_4p_qvec(QList<QVector3D*> coords, double d);
    static arma::mat build_A_3p_qvec(QList<QVector3D*> coords);
    static arma::vec build_b_3p_qvec(QList<QVector3D*> coords, double d);

    static arma::mat build_A_4d(QList<Model::Reservoir::Grid::XYZCoordinate*> coords);
    static arma::vec build_b_4d(QList<Model::Reservoir::Grid::XYZCoordinate*> coords, double d);
    static arma::mat build_A_3d(QList<Model::Reservoir::Grid::XYZCoordinate*> coords);
    static arma::vec build_b_3d(QList<Model::Reservoir::Grid::XYZCoordinate*> coords, double d);
    static arma::vec lambdas(QList<Model::Reservoir::Grid::XYZCoordinate *> coords, double d);
    static Eigen::VectorXd coeff_vector(arma::vec D, arma::mat Qt, arma::vec b);
    static Eigen::Vector3d non_inv_quad_coeffs(Eigen::Vector3d x, Eigen::Vector3d n);

    /* damnit let's try to make everything into Eigen classes. Less errors maybe*/
    static Eigen::Matrix3d build_A_4p_eigen(QList<Eigen::Vector3d> coords);
    static Eigen::Vector3d build_b_4p_eigen(QList<Eigen::Vector3d> coords, double d);
    static Eigen::Matrix3d build_A_3p_eigen(QList<Eigen::Vector3d> coords);
    static Eigen::Vector3d build_b_3p_eigen(QList<Eigen::Vector3d> coords, double d);
    static Eigen::VectorXd coeff_vector_eigen(Eigen::Vector3d D, Eigen::Matrix3d Qt, Eigen::Vector3d b);
    static QList<Eigen::Vector3d> kkt_eq_solutions_eigen(Eigen::Matrix3d A, Eigen::Vector3d b);
    static Eigen::Matrix3d rm_entries_eps_matrix(Eigen::Matrix3d m, double eps);
    static Eigen::VectorXd rm_entries_eps_coeffs(Eigen::VectorXd m, double eps);
    static double shortest_distance_n_wells(QList<QList<Eigen::Vector3d> > coords, int n);
    static double shortest_distance_eigen(QList<Eigen::Vector3d> coords);
    static double shortest_distance_3p_eigen(QList<Eigen::Vector3d> coords);
    static QList<Eigen::Vector3d> well_length_projection_eigen(Eigen::Vector3d heel, Eigen::Vector3d toe, double max, double min, double epsilon);
    static double movement_cost_eig(QList<Eigen::Vector3d> old_coords, QList<Eigen::Vector3d> new_coords);
    static QList<Eigen::Vector3d> move_points_4p_eigen(QList<Eigen::Vector3d> coords, double d, Eigen::Vector3d s);
    static QList<Eigen::Vector3d> move_points_3p_eigen(QList<Eigen::Vector3d> coords, double d, Eigen::Vector3d s);
    static Eigen::Vector3d project_point_to_plane_eigen(Eigen::Vector3d point, Eigen::Vector3d normal_vector, Eigen::Vector3d plane_point);


    /*!
     * \brief Master function for the interwell constraint projection. Projects two wells so that they are at least a distance d appart.
     * \param Two line segments defined by their endpoints. Two first entries belong to first line segment, two last entries belong to second line segment.
     * \return Cheapest(L2-norm) projection of line segments s.t. they are at least a distance d appart.
    */
    static QList<Eigen::Vector3d> interwell_constraint_projection_eigen(QList<Eigen::Vector3d> coords, double d);

    /*!
     * \brief Projects any number of wells so that they are at least a distance d appart.
     * \param Two line segments defined by their endpoints. Two first entries belong to first line segment, two last entries belong to second line segment.
     * \return Cheapest(L2-norm) projection of line segments s.t. they are at least a distance d appart.
    */
    static QList<QList<Eigen::Vector3d> > interwell_constraint_multiple_wells(QList<QList<Eigen::Vector3d> > coords, double d, double tol);
    static void write_eigen_to_tikz(QList<QList<Eigen::Vector3d> > coords, std::string name);
    static QList<QList<Eigen::Vector3d>> well_length_constraint_multiple_wells(QList<QList<Eigen::Vector3d>> wells, double max, double min, double epsilon);
    static bool feasible_well_length(QList<QList<Eigen::Vector3d> > coords, double max, double min, double tol);
    static bool feasible_interwell_distance(QList<QList<Eigen::Vector3d> > coords, double d, double tol);

    static QList<QList<Eigen::Vector3d> > both_constraints_multiple_wells(QList<QList<Eigen::Vector3d> > coords, double d, double tol, double max, double min, double epsilon);


private:
    WellConstraintProjections();


    /*!
     * \brief Sets all elements of a 3-by-3 matrix whose absolute value is below a certain threshold to zero.
     * \param Matrix
     * \param Threshold
     * \return Matrix with cleared sub-threshold values
    */

    static Eigen::Vector3d rm_entries_eps(Eigen::Vector3d m, double eps);

};

}
#endif // WELL_CONSTRAINT_PROJECTIONS_H
