#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/vertex_triangle_adjacency.h>
#include <imgui/imgui.h>
/*** insert any necessary libigl headers here ***/
#include <igl/per_face_normals.h>
#include <igl/copyleft/marching_cubes.h>
#include <math.h>
using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

// Input: imported points, #P x3
Eigen::MatrixXd P;

// Input: imported normals, #P x3
Eigen::MatrixXd N;

// Intermediate result: constrained points, #C x3
Eigen::MatrixXd constrained_points;

// Intermediate result: implicit function values at constrained points, #C x1
Eigen::VectorXd constrained_values;

// Parameter: degree of the polynomial
int deg = 0;

// high value for points with no intersection with constraints
double high = 3.f;

// Parameter: Wendland weight function radius (make this relative to the size of the mesh)
double wendlandRadius = 0.07;

// Parameter: grid resolution
int resolution = 20;

// Intermediate result: grid points, at which the imlicit function will be evaluated, #G x3
Eigen::MatrixXd grid_points;

// Intermediate result: implicit function values at the grid points, #G x1
Eigen::VectorXd grid_values;

// Intermediate result: grid point colors, for display, #G x3
Eigen::MatrixXd grid_colors;

// Intermediate result: grid lines, for display, #L x6 (each row contains
// starting and ending point of line segment)
Eigen::MatrixXd grid_lines;

// Output: vertex array, #V x3
Eigen::MatrixXd V;

// Output: face array, #F x3
Eigen::MatrixXi F;

// Output: face normals of the reconstructed mesh, #F x3
Eigen::MatrixXd FN;

// Expansion ratio to expand grid
double exp_ratio = 20;

// Grid bounds: axis-aligned bounding box
Eigen::RowVector3d bb_min, bb_max;

// Bounding box dimensions
Eigen::RowVector3d dim;

// Grid spacing
double dx, dy, dz;

// Index maps
unordered_map<int, vector<int>> ConstrainedIndex, Pindex;

// Filename
std::string filename;

// PCA alignment
int align_axes = 0;

// Functions
void createGrid();
void getLines();
bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers);
void print_map(const std::unordered_map<int, std::vector<int>> &m);
void nearest(const Eigen::MatrixXd &m, const Eigen::MatrixXd &x, double radius, std::vector<int> &result, const unordered_map<int, vector<int>> &indexMap);
double implicitFunction(const Eigen::MatrixXd &x, double radius, const unordered_map<int, vector<int>> &indexMap);
double implicitFunctionPaper(const Eigen::MatrixXd &x, double radius, const unordered_map<int, vector<int>> &indexMap);
void evaluateImplicitFunc(double radius, const unordered_map<int, vector<int>> &indexMap, std::string version);
void NEAREST(const Eigen::MatrixXd &m, const Eigen::MatrixXd &x, double radius, std::vector<int> &result);

// Creates a grid_points array for the simple sphere example. The points are
// stacked into a single matrix, ordered first in the x, then in the y and
// then in the z direction. If you find it necessary, replace this with your own
// function for creating the grid.
void createGrid()
{
    grid_points.resize(0, 3);
    grid_colors.resize(0, 3);
    grid_lines.resize(0, 6);
    grid_values.resize(0);
    V.resize(0, 3);
    F.resize(0, 3);
    FN.resize(0, 3);

    // 3D positions of the grid points -- see slides or marching_cubes.h for ordering
    grid_points.resize(resolution * resolution * resolution, 3);
    // Create each gridpoint
    for (unsigned int x = 0; x < resolution; ++x)
    {
        for (unsigned int y = 0; y < resolution; ++y)
        {
            for (unsigned int z = 0; z < resolution; ++z)
            {
                // Linear index of the point at (x,y,z)
                int index = x + resolution * (y + resolution * z);
                // 3D point at (x,y,z)
                grid_points.row(index) = bb_min + Eigen::RowVector3d(x * dx, y * dy, z * dz);
            }
        }
    }
}

// Function for explicitly evaluating the implicit function for a sphere of
// radius r centered at c : f(p) = ||p-c|| - r, where p = (x,y,z).
// This will NOT produce valid results for any mesh other than the given
// sphere.
// Replace this with your own function for evaluating the implicit function
// values at the grid points using MLS
void evaluateImplicitFuncSphere()
{
    // Sphere center
    auto bb_min = grid_points.colwise().minCoeff().eval();
    auto bb_max = grid_points.colwise().maxCoeff().eval();
    Eigen::RowVector3d center = 0.5 * (bb_min + bb_max);

    double radius = 0.5 * (bb_max - bb_min).minCoeff();

    // Scalar values of the grid points (the implicit function values)
    grid_values.resize(resolution * resolution * resolution);

    // Evaluate sphere's signed distance function at each gridpoint.
    for (unsigned int x = 0; x < resolution; ++x)
    {
        for (unsigned int y = 0; y < resolution; ++y)
        {
            for (unsigned int z = 0; z < resolution; ++z)
            {
                // Linear index of the point at (x,y,z)
                int index = x + resolution * (y + resolution * z);

                // Value at (x,y,z) = implicit function for the sphere
                grid_values[index] = (grid_points.row(index) - center).norm() - radius;
            }
        }
    }
}

void evaluateImplicitFunc(double radius, const std::unordered_map<int, std::vector<int>> &indexMap, std::string version)
{
    // Scalar values of the grid points (the implicit function values)
    grid_values.resize(resolution * resolution * resolution);

    // Evaluate signed distance function at each gridpoint.
    for (unsigned int x = 0; x < resolution; ++x)
    {
        for (unsigned int y = 0; y < resolution; ++y)
        {
            for (unsigned int z = 0; z < resolution; ++z)
            {
                // Linear index of the point at (x,y,z)
                int index = x + resolution * (y + resolution * z);
                Eigen::MatrixXd X(1, 3);
                X.row(0) << grid_points.row(index);
                if (version == "standard")
                {
                    grid_values[index] = implicitFunction(X, radius, indexMap);
                }
                else if (version == "paper")
                {
                    grid_values[index] = implicitFunctionPaper(X, radius, indexMap);
                }
                else
                {
                    cout << "Version not recognized! Use \"standard\" or \"paper\"" << endl;
                }
            }
        }
    }
}

Eigen::MatrixXd poly_vector(const Eigen::MatrixXd &p)
{
    // Assume the degree can be 0, 1 or 2
    Eigen::MatrixXd res;
    const int n = p.rows(); // Number of constrained points
    if (deg == 0)
    {
        res.resize(n, 1);
        res.setOnes(n, 1);
        return res;
    }
    else if (deg == 1)
    {
        res.resize(n, 4);
        res.setOnes(n, 4);
        res.block(0, 1, p.rows(), 3) = p;
        return res;
    }
    else
    {
        res.resize(n, 10);
        res.setOnes(n, 10);
        res.block(0, 1, p.rows(), p.cols()) = p;
        res.block(0, 4, p.rows(), p.cols()) = p.cwiseAbs2();
        res.block(0, 7, p.rows(), 1) = p.col(0).cwiseProduct(p.col(1));
        res.block(0, 8, p.rows(), 1) = p.col(1).cwiseProduct(p.col(2));
        res.block(0, 9, p.rows(), 1) = p.col(2).cwiseProduct(p.col(0));
        return res;
    }
}

Eigen::MatrixXd proximity_weights_matrix(const Eigen::MatrixXd &x, const Eigen::MatrixXd &p, const double radius)
{
    Eigen::MatrixXd v;
    v = p;
    v.rowwise() -= Eigen::RowVectorXd(x);
    Eigen::MatrixXd val(p.rows(), 1);
    val = v.rowwise().norm();
    // Apply wendland on val
    val /= radius;
    val = ((1 - val.array()).pow(4) * ((4 * val).array() + 1)).matrix();
    val = Eigen::VectorXd(val).asDiagonal();
    return val.cwiseSqrt();
}

double implicitFunction(const Eigen::MatrixXd &x, double radius, const std::unordered_map<int, std::vector<int>> &indexMap)
{
    // high is an arbitrary high value assigned to points in whose neighborhood there are
    // less than the min number of points to fit the polynomial
    // Find neighbors of x within radius
    std::vector<int> neighbors;
    //NEAREST(constrained_points, x, radius, neighbors);
    nearest(constrained_points, x, radius, neighbors, indexMap);
    std::unordered_map<int, int> dict;
    dict[0] = 1;
    dict[1] = 4;  // 4
    dict[2] = 10; // 10
    if (neighbors.size() < dict[deg])
    {
        return high;
    }
    else
    {
        // Build p
        Eigen::MatrixXd p(neighbors.size(), 3);
        for (int i = 0; i < neighbors.size(); i++)
        {
            p.row(i) << constrained_points.row(neighbors[i]);
        }
        Eigen::MatrixXd proximity_w;
        proximity_w = proximity_weights_matrix(x, p, radius);
        Eigen::MatrixXd poly_vect;
        poly_vect = poly_vector(p);
        // Extract constrained values
        Eigen::VectorXd constraints(neighbors.size(), 1);
        for (int i = 0; i < neighbors.size(); i++)
        {
            constraints(i) = constrained_values(neighbors[i]);
        }
        // Find the coefficients
        Eigen::MatrixXd A;
        Eigen::VectorXd b;
        A = proximity_w * poly_vect;
        b = proximity_w * constraints;
        Eigen::VectorXd c = A.colPivHouseholderQr().solve(b);
        //Eigen::VectorXd c = (A.transpose() * A).ldlt().solve(A.transpose() * b);

        // Compute solution
        double computed_result;
        Eigen::MatrixXd basis_x;
        basis_x = poly_vector(x);
        Eigen::VectorXd vectorized_basis;
        vectorized_basis = basis_x.row(0);
        computed_result = (vectorized_basis).dot(c);
        return computed_result;
    }
}

double implicitFunctionPaper(const Eigen::MatrixXd &x, double radius, const std::unordered_map<int, std::vector<int>> &indexMap)
{
    // high is an arbitrary high value assigned to points in whose neighborhood there are
    // less than the min number of points to fit the polynomial
    // Find neighbors of x within radius
    std::vector<int> neighbors;
    nearest(P, x, radius, neighbors, indexMap);
    std::unordered_map<int, int> dict;
    dict[0] = 1;
    dict[1] = 4; // 4
    dict[2] = 10;
    if (neighbors.size() < dict[deg])
    {
        return high;
    }
    else
    {
        // Build p
        Eigen::MatrixXd p(neighbors.size(), 3);
        for (int i = 0; i < neighbors.size(); i++)
        {
            p.row(i) << P.row(neighbors[i]);
        }
        Eigen::MatrixXd proximity_w;
        proximity_w = proximity_weights_matrix(x, p, radius);
        Eigen::MatrixXd poly_vect;
        poly_vect = poly_vector(p);
        // Compute constrained values
        Eigen::VectorXd constraints(neighbors.size(), 1);
        Eigen::MatrixXd N_norm(P.rows(), 3);
        //N_norm.setZero(P.rows() * 3, 3);
        N_norm = N.rowwise().normalized();

        for (int i = 0; i < neighbors.size(); i++)
        {
            constraints(i) = ((Eigen::Vector3d)(x - P.row(neighbors[i])).transpose()).dot((Eigen::Vector3d)N_norm.row(neighbors[i]).transpose());
        }
        // Find the coefficients
        Eigen::MatrixXd A;
        Eigen::VectorXd b;
        A = proximity_w * poly_vect;
        b = proximity_w * constraints;
        Eigen::VectorXd c = A.colPivHouseholderQr().solve(b);
        //Eigen::VectorXd c = (A.transpose() * A).ldlt().solve(A.transpose() * b);
        // Compute solution
        double computed_result;
        Eigen::MatrixXd basis_x;
        basis_x = poly_vector(x);
        Eigen::VectorXd vectorized_basis;
        vectorized_basis = basis_x.row(0);
        computed_result = (vectorized_basis).dot(c);
        return computed_result;
    }
}

// Code to display the grid lines given a grid structure of the given form.
// Assumes grid_points have been correctly assigned
// Replace with your own code for displaying lines if need be.
void getLines()
{
    int nnodes = grid_points.rows();
    grid_lines.resize(3 * nnodes, 6);
    int numLines = 0;

    for (unsigned int x = 0; x < resolution; ++x)
    {
        for (unsigned int y = 0; y < resolution; ++y)
        {
            for (unsigned int z = 0; z < resolution; ++z)
            {
                int index = x + resolution * (y + resolution * z);
                if (x < resolution - 1)
                {
                    int index1 = (x + 1) + y * resolution + z * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (y < resolution - 1)
                {
                    int index1 = x + (y + 1) * resolution + z * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (z < resolution - 1)
                {
                    int index1 = x + y * resolution + (z + 1) * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
            }
        }
    }

    grid_lines.conservativeResize(numLines, Eigen::NoChange);
}

int get_box(Eigen::MatrixXd x)
{
    // Given a 3d point, return the index of the grid box it is in
    x = x - bb_min;
    x << min((int)floor(x(0) / dx), resolution - 2), min((int)floor(x(1) / dy), resolution - 2), min((int)floor(x(2) / dz), resolution - 2);
    int ind = x(0, 0) + resolution * (x(0, 1) + resolution * x(0, 2));
    return ind;
}

Eigen::RowVector3i get_integer_coordinates(Eigen::MatrixXd x)
{
    // Given a 3d point, return the integer coordinates of the grid box representative vertex it is in
    x = x - bb_min;
    Eigen::RowVector3i res;
    res << min((int)floor(x(0, 0) / dx), resolution - 2), min((int)floor(x(0, 1) / dy), resolution - 2), min((int)floor(x(0, 2) / dz), resolution - 2);
    return res;
}

void buildIndex(Eigen::MatrixXd &m, std::unordered_map<int, std::vector<int>> &indexMap)
{
    // We build a dictionary where the key is the index of a vertex of the grid (a vertex
    // represents the cube that has that vertex as lower, left, ner to observer), and the
    // value is a list of all indices of points P in that cube
    indexMap.clear();
    for (int r = 0; r < m.rows(); ++r)
    {
        // Find vertex that indentifies cube in which m.row lies
        int ind = get_box(m.row(r));
        auto it = indexMap.find(ind);
        if (it == indexMap.end())
        {
            indexMap[ind] = {r};
        }
        else
        {
            it->second.push_back(r);
        }
    }
}

int integer_coord_to_index(Eigen::RowVector3i x)
{
    int ind = x(0, 0) + resolution * (x(0, 1) + resolution * x(0, 2));
    return ind;
}

void nearest(const Eigen::MatrixXd &m, const Eigen::MatrixXd &x, double radius, std::vector<int> &result, const unordered_map<int, std::vector<int>> &indexMap)
{
    // Return all points of m within distance "radius" from x, ordered from nearer to farthest
    std::vector<pair<double, int>> intermediate;
    // Empty result
    result.clear();
    Eigen::RowVector3i point, new_location;
    point = get_integer_coordinates(x);

    for (int i_x = -ceil(radius / dx) - 1; i_x < ceil(radius / dx) + 1; ++i_x)
    {
        for (int i_y = -ceil(radius / dy) - 1; i_y < ceil(radius / dy) + 1; ++i_y)
        {
            for (int i_z = -ceil(radius / dz) - 1; i_z < ceil(radius / dz) + 1; ++i_z)
            {
                Eigen::RowVector3i direction = {i_x, i_y, i_z};
                new_location = point + direction;
                if ((new_location.array() >= 0).all() && (new_location.array() < resolution - 1).all())
                {
                    int ind = integer_coord_to_index(new_location);
                    const auto it = indexMap.find(ind);
                    if (it != indexMap.end())
                    {
                        for (const auto el : it->second)
                        {
                            double dist = (m.row(el) - x).norm();
                            if (dist <= radius)
                            {
                                intermediate.push_back(make_pair(dist, el));
                            }
                        }
                    }
                }
            }
        }
    }
    std::sort(intermediate.begin(), intermediate.end());
    // Build final result
    for (int i = 0; i < intermediate.size(); i++)
    {
        result.push_back(intermediate[i].second);
    }
}

void find_epsilon(double init_epsilon, double &epsilon, Eigen::MatrixXd N_norm, const std::unordered_map<int, std::vector<int>> &indexMap)
{
    cout << "start find eps" << endl;
    std::vector<int> near_out = {}, near_in = {};

    for (int i = 0; i < P.rows(); ++i)
    {
        nearest(P, P.row(i) + (N_norm.row(i) * init_epsilon), init_epsilon * 3, near_out, indexMap);
        nearest(P, P.row(i) - (N_norm.row(i) * init_epsilon), init_epsilon * 3, near_in, indexMap);
        while (near_out[0] != i || near_in[0] != i)
        {
            init_epsilon = init_epsilon / 2;
            nearest(P, P.row(i) + (N_norm.row(i) * init_epsilon), init_epsilon * 3, near_out, indexMap);
            nearest(P, P.row(i) - (N_norm.row(i) * init_epsilon), init_epsilon * 3, near_in, indexMap);
        }
    }
    epsilon = init_epsilon;
}

void NEAREST(const Eigen::MatrixXd &m, const Eigen::MatrixXd &x, double radius, std::vector<int> &result)
{
    // Clear
    result.clear();
    std::vector<pair<double, int>> intermediate;
    double dist;
    for (int i = 0; i < m.rows(); i++)
    {
        dist = (m.row(i) - x).norm();
        if (dist < radius)
        {
            intermediate.push_back(make_pair(dist, i));
        }
    }
    std::sort(intermediate.begin(), intermediate.end());
    // Build final result
    for (int i = 0; i < intermediate.size(); i++)
    {
        result.push_back(intermediate[i].second);
    }
}

void align_pointcloud(Eigen::MatrixXd &mat, Eigen::MatrixXd &normals)
{
    Eigen::MatrixXd centered = mat.rowwise() - mat.colwise().mean();
    Eigen::MatrixXd cov = centered.adjoint() * centered;
    // Eigen decomposition
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);
    Eigen::Vector3d most_variant_axis;
    most_variant_axis = eig.eigenvectors().rightCols(1);
    // find rotation to align direction with highest variance to axis y
    Eigen::Vector3d axis;
    axis = most_variant_axis.cross(Eigen::Vector3d::UnitY());
    double angle;
    angle = asin(axis.norm() / (most_variant_axis.norm()));
    Eigen::AngleAxisd angleAxis(angle, axis.normalized());
    Eigen::Matrix3d rotation = angleAxis.toRotationMatrix();
    // Perform rotation of pointcloud (modify it in place)
    mat = (rotation * mat.transpose()).transpose();
    normals = (rotation * normals.transpose()).transpose();
}

void print_map(const std::unordered_map<int, std::vector<int>> &m)
{
    for (const auto &pair : m)
    {
        std::cout << "{" << pair.first << ": [";
        for (auto el : pair.second)
        {
            cout << el << ", ";
        }
        cout << "] }" << endl;
    }
}

void print_vect(const std::vector<int> &v)
{
    std::cout << "[ ";
    for (auto el : v)
    {
        std::cout << el << ", ";
    }
    std::cout << "]" << endl;
}

void debug_nearest(const Eigen::MatrixXd &x, const Eigen::MatrixXd &m, double radius, const std::unordered_map<int, std::vector<int>> &indexMap)
{
    std::vector<int> result1 = {}, result2 = {};
    std::cout << "Ground truth: " << endl;
    NEAREST(m, x, radius, result2);
    print_vect(result2);
    std::cout << "nearest result: " << endl;
    nearest(m, x, radius, result1, indexMap);
    print_vect(result1);
}

void build_constraints(const Eigen::MatrixXd &N, const unordered_map<int, vector<int>> &indexMap)
{
    // Normalized normals
    Eigen::MatrixXd N_norm(P.rows(), 3);
    N_norm = N.rowwise().normalized();

    // n: number of vertices
    int n = P.rows();

    // Resize
    constrained_points.resize(3 * n, 3);
    constrained_values.resize(3 * n);

    // Clear
    constrained_points.setZero(3 * n, 3);
    constrained_values.setZero(3 * n);

    // Add first n simpler points and constraints
    constrained_points.block(0, 0, n, 3) = P;

    // Find appropriate epsilon
    double diagonal = dim.norm();
    double init_epsilon = 0.03 * diagonal;
    double epsilon;

    find_epsilon(init_epsilon, epsilon, N_norm, indexMap);

    // Add outer/inner off surface constraints
    constrained_points.block(n, 0, n, 3) = P + (N_norm * epsilon);
    constrained_points.block(2 * n, 0, n, 3) = P - (N_norm * epsilon);
    for (int i = 0; i < n; i++)
    {
        constrained_values(n + i) = epsilon;
        constrained_values(2 * n + i) = -epsilon;
    }
}

bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers)
{
    if (key == '1')
    {
        // Show imported points
        viewer.data().clear();
        viewer.core.align_camera_center(P);
        viewer.data().point_size = 11;
        viewer.data().add_points(P, Eigen::RowVector3d(0, 0, 0));
        // Test
        /*igl::readOFF(filename, V, F);
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        viewer.data().compute_normals();
        viewer.core.align_camera_center(V, F);
        viewer.data().clear();
        viewer.data().set_mesh(V, F);*/
    }

    if (key == '2')
    {
        // Build color map
        Eigen::MatrixXd C;
        C.setZero(P.rows() * 3, 3);
        for (int i = 0; i < C.rows(); ++i)
        {

            if (i < P.rows())
            {
                C(i, 0) = 1;
            }
            else if (i < 2 * P.rows())
            {

                C(i, 1) = 1;
            }
            else
            {
                C(i, 2) = 1;
            }
        }
        // Show imported points
        viewer.data().clear();
        viewer.core.align_camera_center(P);
        viewer.data().point_size = 6;
        viewer.data().add_points(constrained_points, C);
    }

    if (key == '3')
    {

        // Show grid points with colored nodes and connected with lines
        viewer.data().clear();
        viewer.core.align_camera_center(P);
        // Add code for creating a grid
        // Add your code for evaluating the implicit function at the grid points
        // Add code for displaying points and lines
        // You can use the following example:

        /*** begin: sphere example, replace (at least partially) with your code ***/
        // Make grid
        createGrid();

        // Evaluate implicit function
        evaluateImplicitFunc(wendlandRadius * dim.norm(), ConstrainedIndex, "standard");

        // get grid lines
        getLines();

        // Code for coloring and displaying the grid points and lines
        // Assumes that grid_values and grid_points have been correctly assigned.
        grid_colors.setZero(grid_points.rows(), 3);

        // Build color map
        for (int i = 0; i < grid_points.rows(); ++i)
        {
            double value = grid_values(i);
            if (value < 0)
            {
                grid_colors(i, 1) = 1;
            }
            else
            {
                if (value > 0)
                    grid_colors(i, 0) = 1;
            }
        }

        // Draw lines and points
        viewer.data().point_size = 8;
        viewer.data().add_points(grid_points, grid_colors);
        viewer.data().add_edges(grid_lines.block(0, 0, grid_lines.rows(), 3),
                                grid_lines.block(0, 3, grid_lines.rows(), 3),
                                Eigen::RowVector3d(0.8, 0.8, 0.8));
        /*** end: sphere example ***/
    }

    if (key == '4')
    {
        // Make grid
        createGrid();

        // Evaluate implicit function
        evaluateImplicitFunc(wendlandRadius * dim.norm(), ConstrainedIndex, "standard");

        // Show reconstructed mesh
        viewer.data().clear();
        // Code for computing the mesh (V,F) from grid_points and grid_values
        if ((grid_points.rows() == 0) || (grid_values.rows() == 0))
        {
            cerr << "Not enough data for Marching Cubes !" << endl;
            return true;
        }
        // Run marching cubes
        igl::copyleft::marching_cubes(grid_values, grid_points, resolution, resolution, resolution, V, F);
        if (V.rows() == 0)
        {
            cerr << "Marching Cubes failed!" << endl;
            return true;
        }

        igl::per_face_normals(V, F, FN);
        viewer.data().set_mesh(V, F);
        viewer.data().show_lines = true;
        viewer.data().show_faces = true;
        viewer.data().set_normals(FN);

        // Write mesh
        std::string s = "../results/reconstructed_" + filename.substr(8, filename.length() - 7);
        igl::writeOFF(s, V, F);
    }

    if (key == '5')
    {
        // Make grid
        createGrid();

        // Evaluate implicit function
        evaluateImplicitFunc(wendlandRadius * dim.norm(), Pindex, "paper");

        // Show reconstructed mesh
        viewer.data().clear();
        // Code for computing the mesh (V,F) from grid_points and grid_values
        if ((grid_points.rows() == 0) || (grid_values.rows() == 0))
        {
            cerr << "Not enough data for Marching Cubes !" << endl;
            return true;
        }
        // Run marching cubes
        igl::copyleft::marching_cubes(grid_values, grid_points, resolution, resolution, resolution, V, F);
        if (V.rows() == 0)
        {
            cerr << "Marching Cubes failed!" << endl;
            return true;
        }

        igl::per_face_normals(V, F, FN);
        viewer.data().set_mesh(V, F);
        viewer.data().show_lines = true;
        viewer.data().show_faces = true;
        viewer.data().set_normals(FN);
    }

    return true;
}

bool callback_load_mesh(Viewer &viewer, string filename)
{
    igl::readOFF(filename, P, F, N);
    callback_key_down(viewer, '1', 0);
    return true;
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cout << "Usage ex2_bin <mesh.off>" << endl;
        filename = "../data/sphere.off";
        igl::readOFF(filename, P, F, N);
    }
    else
    {
        // Read points and normals
        igl::readOFF(argv[1], P, F, N);
        filename = argv[1];
    }

    // Align the pointcloud with the axes
    if (filename == "../data/luigi.off")
    {
        align_axes = 1;
        align_pointcloud(P, N);
    }

    if (filename == "../data/luigi.off" || filename == "../data/hound.off" || filename == "../data/horse.off")
    {
        resolution = 40;
        wendlandRadius = 0.02;
    }

    // Grid bounds: axis-aligned bounding box
    bb_min = P.colwise().minCoeff();
    bb_max = P.colwise().maxCoeff();

    // Bounding box dimensions
    dim = bb_max - bb_min;

    // Enlarge the grid!
    bb_min = bb_min.array() - dim.norm() / exp_ratio;
    bb_max = bb_max.array() + dim.norm() / exp_ratio;
    dim = bb_max - bb_min;

    // Grid spacing
    dx = dim[0] / (double)(resolution - 1);
    dy = dim[1] / (double)(resolution - 1);
    dz = dim[2] / (double)(resolution - 1);

    // Build the constrained points and relative constraints
    buildIndex(P, Pindex);
    build_constraints(N, Pindex);

    // Map for constraints
    buildIndex(constrained_points, ConstrainedIndex);

    Viewer viewer;
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

    viewer.callback_key_down = callback_key_down;

    menu.callback_draw_viewer_menu = [&]() {
        // Draw parent menu content
        menu.draw_viewer_menu();

        // Add new group
        if (ImGui::CollapsingHeader("Reconstruction Options", ImGuiTreeNodeFlags_DefaultOpen))
        {

            // Degree of polynomial
            ImGui::InputInt("Poly deg", &deg, 0, 0);
            if (ImGui::Button("Poly degree", ImVec2(-1, 0)))
            {
                cout << "Poly degree set to " << deg << endl;
                // Switch view to show paper implementation with new radius
                //callback_key_down(viewer, '4', 0);
            }

            // Resolution
            ImGui::InputInt("Resolution", &resolution, 0, 0);
            if (ImGui::Button("Reset Grid", ImVec2(-1, 0)))
            {
                std::cout << "Resolution set to " << resolution << endl;

                // Recreate the grid
                // createGrid();
                // Grid spacing
                dx = dim[0] / (double)(resolution - 1);
                dy = dim[1] / (double)(resolution - 1);
                dz = dim[2] / (double)(resolution - 1);

                // Build the constrained points and relative constraints
                buildIndex(P, Pindex);
                build_constraints(N, Pindex);

                // Map for constraints
                buildIndex(constrained_points, ConstrainedIndex);
            }
            // Button to align exes
            if (ImGui::Button("Align axes", ImVec2(-1, 0)))
            {
                std::cout << "Axes aligned " << endl;

                align_pointcloud(P, N);

                // Grid bounds: axis-aligned bounding box
                bb_min = P.colwise().minCoeff();
                bb_max = P.colwise().maxCoeff();

                // Bounding box dimensions
                dim = bb_max - bb_min;

                // Enlarge the grid!
                bb_min = bb_min.array() - dim.norm() / exp_ratio;
                bb_max = bb_max.array() + dim.norm() / exp_ratio;
                dim = bb_max - bb_min;

                // Grid spacing
                dx = dim[0] / (double)(resolution - 1);
                dy = dim[1] / (double)(resolution - 1);
                dz = dim[2] / (double)(resolution - 1);

                // Build the constrained points and relative constraints
                buildIndex(P, Pindex);
                build_constraints(N, Pindex);

                // Map for constraints
                buildIndex(constrained_points, ConstrainedIndex);
                callback_key_down(viewer, '1', 0);
            }

            // Wendland radius
            ImGui::InputDouble("Wendland Radius", &wendlandRadius, 0, 0);
            if (ImGui::Button("Modify Wendland Radius", ImVec2(-1, 0)))
            {
                cout << "Wendland radius set to " << wendlandRadius << endl;
                // Switch view to show paper implementation with new radius
                //callback_key_down(viewer, '4', 0);
            }
        }
    };

    callback_key_down(viewer, '1', 0);

    viewer.launch();
}
