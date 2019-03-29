#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
/*** insert any necessary libigl headers here ***/
#include <igl/per_face_normals.h>
#include <igl/copyleft/marching_cubes.h>

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
int polyDegree = 0;

// Parameter: Wendland weight function radius (make this relative to the size of the mesh)
double wendlandRadius = 0.1;

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

// Functions
void createGrid();
void evaluateImplicitFunc();
void getLines();
bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers);

// Creates a grid_points array for the simple sphere example. The points are
// stacked into a single matrix, ordered first in the x, then in the y and
// then in the z direction. If you find it necessary, replace this with your own
// function for creating the grid.
void createGrid() {
    grid_points.resize(0, 3);
    grid_colors.resize(0, 3);
    grid_lines. resize(0, 6);
    grid_values.resize(0);
    V. resize(0, 3);
    F. resize(0, 3);
    FN.resize(0, 3);

    // Grid bounds: axis-aligned bounding box
    Eigen::RowVector3d bb_min, bb_max;
    bb_min = P.colwise().minCoeff();
    bb_max = P.colwise().maxCoeff();

    // Bounding box dimensions
    Eigen::RowVector3d dim = bb_max - bb_min;

    // Grid spacing
    const double dx = dim[0] / (double)(resolution - 1);
    const double dy = dim[1] / (double)(resolution - 1);
    const double dz = dim[2] / (double)(resolution - 1);

    // 3D positions of the grid points -- see slides or marching_cubes.h for ordering
    grid_points.resize(resolution * resolution * resolution, 3);
    // Create each gridpoint
    for (unsigned int x = 0; x < resolution; ++x) {
        for (unsigned int y = 0; y < resolution; ++y) {
            for (unsigned int z = 0; z < resolution; ++z) {
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
void evaluateImplicitFunc() {

    // Scalar values of the grid points (the implicit function values)

    grid_values.resize(resolution * resolution * resolution);
    
    // Radius of neighbors

    double radius = 4;

    // Evaluate sphere's signed distance function at each gridpoint.

    for (unsigned int x = 0; x < resolution; ++x) {
        for (unsigned int y = 0; y < resolution; ++y) {
            for (unsigned int z = 0; z < resolution; ++z) {

                // Linear index of the point at (x,y,z)

                int index = x + resolution * (y + resolution * z);

                // Current grid point

                Eigen::RowVector3d point1;
                point1(0) = grid_points(index, 0);
                point1(1) = grid_points(index, 1);
                point1(2) = grid_points(index, 2);

                // Initialization

                Eigen::MatrixXd A;
                Eigen::MatrixXd f;
                Eigen::MatrixXd w;

                A.resize(0, 4);
                f.resize(0, 1);
                w.resize(0, 0);

                // Searching closest neighbors

                for (int i = 0; i < constrained_points.rows(); i++) {

                    Eigen::RowVector3d point3;

                    point3(0) = constrained_points(i, 0);
                    point3(1) = constrained_points(i, 1);
                    point3(2) = constrained_points(i, 2);

                    Eigen::RowVector3d diff = point1 - point3;

                    double distance = sqrt(pow(diff(0), 2) + pow(diff(1), 2) + pow(diff(2), 2));

                    if (distance < radius) {

                        A.conservativeResize(A.rows() + 1, 4);
                        f.conservativeResize(f.rows() + 1, 1);
                        w.conservativeResize(w.rows() + 1, w.rows() + 1);

                        // Update A

                        A(A.rows()-1, 0) = 1;
                        for (int j = 0; j < 3; j++) {
                            A(A.rows() - 1, j + 1) = constrained_points(i, j);
                        }

                        // Update f

                        f(f.rows()-1, 0) = constrained_values(i, 0);

                        // Update w

                        for (int j = 0; j < w.rows(); j++) {
                            w(j, w.rows() - 1) = 0;
                            w(w.rows() - 1, j) = 0;
                        }

                        w(w.rows()-1, w.rows() - 1) = distance;

                    }

                }

                // Solving the linear equations

                Eigen::MatrixXd A_tag = w * A;
                Eigen::MatrixXd f_tag = w * f;
                Eigen::RowVector4d c = A.colPivHouseholderQr().solve(f_tag).transpose();

                // value of current grid point

                Eigen::RowVector4d point2;
                point2(0) = 1;
                point2(1) = grid_points(index, 0);
                point2(2) = grid_points(index, 1);
                point2(3) = grid_points(index, 2);

                grid_values[index] = point2.dot(c);

            }
        }
    }

}

// Code to display the grid lines given a grid structure of the given form.
// Assumes grid_points have been correctly assigned
// Replace with your own code for displaying lines if need be.
void getLines() {
    int nnodes = grid_points.rows();
    grid_lines.resize(3 * nnodes, 6);
    int numLines = 0;

    for (unsigned int x = 0; x<resolution; ++x) {
        for (unsigned int y = 0; y < resolution; ++y) {
            for (unsigned int z = 0; z < resolution; ++z) {
                int index = x + resolution * (y + resolution * z);
                if (x < resolution - 1) {
                    int index1 = (x + 1) + y * resolution + z * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (y < resolution - 1) {
                    int index1 = x + (y + 1) * resolution + z * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (z < resolution - 1) {
                    int index1 = x + y * resolution + (z + 1) * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
            }
        }
    }

    grid_lines.conservativeResize(numLines, Eigen::NoChange);
}

bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers) {
    if (key == '1') {
        // Show imported points
        viewer.data().clear();
        viewer.core.align_camera_center(P);

        viewer.data().point_size = 11;
        viewer.data().add_points(P, Eigen::RowVector3d(0,0,0));
    }

    if (key == '2') {
        // Show all constraints
        viewer.data().clear();
        viewer.core.align_camera_center(P);

        // Add your code for computing auxiliary constraint points here

        // Bounding box dimensions
        Eigen::RowVector3d bb_min, bb_max;
        bb_min = P.colwise().minCoeff();
        bb_max = P.colwise().maxCoeff();
        Eigen::RowVector3d dim = bb_max - bb_min;

        double epsilon = 0.01 * sqrt(pow(dim[0], 2) + pow(dim[1], 2) + pow(dim[2], 2));

        constrained_points.resize(3 * P.rows(), 3);
        constrained_values.resize(3 * P.rows(), 1);

        Eigen::MatrixXd colors_per_face;
        colors_per_face.setZero(3 * P.rows(), 3);

        for (int i = 0; i < P.rows(); i++) {

            constrained_points.row(3 * i)       = P.row(i);
            constrained_values(3 * i, 0)        = 0;
            colors_per_face.row(3 * i)          = Eigen::RowVector3d(0, 0, 1);

            constrained_points.row(3 * i + 1)   = P.row(i) + epsilon * N.row(i);
            constrained_values(3 * i + 1, 0)    = epsilon;
            colors_per_face.row(3 * i + 1)      = Eigen::RowVector3d(1, 0, 0);

            constrained_points.row(3 * i + 2)   = P.row(i) - epsilon * N.row(i);
            constrained_values(3 * i + 2, 0)    = -epsilon;
            colors_per_face.row(3 * i + 2)      = Eigen::RowVector3d(0, 1, 0);

        }

        // Add code for displaying all points, as above

        viewer.data().point_size = 11;
        viewer.data().add_points(constrained_points, colors_per_face);

    }

    if (key == '3') {
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
        evaluateImplicitFunc();

        // get grid lines
//        getLines();

        // Code for coloring and displaying the grid points and lines
        // Assumes that grid_values and grid_points have been correctly assigned.
        grid_colors.setZero(grid_points.rows(), 3);

        // Build color map
        for (int i = 0; i < grid_points.rows(); ++i) {
            double value = grid_values(i);
            if (value < 0) {
                grid_colors(i, 1) = 1;
            }
            else {
                if (value > 0)
                    grid_colors(i, 0) = 1;
            }
        }

        // Draw lines and points
        viewer.data().point_size = 8;
        viewer.data().add_points(grid_points, grid_colors);

//       viewer.data().add_edges(grid_lines.block(0, 0, grid_lines.rows(), 3), grid_lines.block(0, 3, grid_lines.rows(), 3), Eigen::RowVector3d(0.8, 0.8, 0.8));

        /*** end: sphere example ***/
    }

    if (key == '4') {
        // Show reconstructed mesh
        viewer.data().clear();
        // Code for computing the mesh (V,F) from grid_points and grid_values
        if ((grid_points.rows() == 0) || (grid_values.rows() == 0)) {
            cerr << "Not enough data for Marching Cubes !" << endl;
            return true;
        }
        // Run marching cubes
        igl::copyleft::marching_cubes(grid_values, grid_points, resolution, resolution, resolution, V, F);
        if (V.rows() == 0) {
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

bool callback_load_mesh(Viewer& viewer,string filename)
{
  igl::readOFF(filename,P,F,N);
  callback_key_down(viewer,'1',0);
  return true;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
      cout << "Usage ex2_bin <mesh.off>" << endl;
      igl::readOFF("../data/Test.off", P, F, N);
//      igl::readOFF("../data/sphere.off",P,F,N);
//      igl::readOFF("../data/cat.off", P, F, N);
//      igl::readOFF("../data/bunny-500.off", P, F, N);
    }
	  else
	  {
		  // Read points and normals
		  igl::readOFF(argv[1],P,F,N);
	  }

    Viewer viewer;
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

    viewer.callback_key_down = callback_key_down;

    menu.callback_draw_viewer_menu = [&]()
    {
      // Draw parent menu content
      menu.draw_viewer_menu();

      // Add new group
      if (ImGui::CollapsingHeader("Reconstruction Options", ImGuiTreeNodeFlags_DefaultOpen))
      {
        // Expose variable directly ...
        ImGui::InputInt("Resolution", &resolution, 0, 0);
        if (ImGui::Button("Reset Grid", ImVec2(-1,0)))
        {
          std::cout << "ResetGrid\n";
          // Recreate the grid
          createGrid();
          // Switch view to show the grid
          callback_key_down(viewer,'3',0);
        }

        // TODO: Add more parameters to tweak here...
      }

    };

    callback_key_down(viewer, '1', 0);

    viewer.launch();
}
