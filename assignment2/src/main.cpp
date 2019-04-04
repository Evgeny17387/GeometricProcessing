#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
/*** insert any necessary libigl headers here ***/
#include <igl/per_face_normals.h>
#include <igl/copyleft/marching_cubes.h>

#define MAX_INT 0xFFFFFFFF
#define EPSILON 0.00000001

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
int polyDegree = 1;

// Parameter: Wendland weight function radius (make this relative to the size of the mesh)
double wendlandRadius_factor = 0.1;

// Parameter: Factor multiplies object dimensions for epsilon
double epsilon_factor = 0.01;

// Parameter: Smoothing factor for weight function
double h = 1;

// Parameter: Object number
int object = 0;

// Parameter: grid resolution
int resolutionX = 10;
int resolutionY = 10;
int resolutionZ = 10;

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

// another grid parameters
list<int> * constrained_grid_points = NULL;
int resolutionXRaw = 10;
int resolutionYRaw = 10;
int resolutionZRaw = 10;

// rotation angles


// Functions
void rotate();
double weightFunction(double r);
void createGrid();
void evaluateImplicitFunc();
void getLines();
bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers);

double theta = 0;
double phi = 0;
double zeta = 0;

#define PI 3.14

// Rotates Vertices and Normals
void rotate() {

    double cTheta = cos(theta*PI / 180);
    double sTheta = sin(theta*PI / 180);

    Eigen::Matrix3d Rx;

    Rx  <<  1, 0, 0,
            0, cTheta, -sTheta,
            0, sTheta, cTheta;

    double cPhi = cos(phi*PI / 180);
    double sPhi = sin(phi*PI / 180);

    Eigen::Matrix3d Ry;

    Ry  <<  cPhi, 0, sPhi,
            0, 1, 0,
            -sPhi, 0, cPhi;

    double cZeta = cos(zeta*PI / 180);
    double sZeta = sin(zeta*PI / 180);

    Eigen::Matrix3d Rz;

    Rz  <<  cZeta, -sZeta, 0,
            sZeta, cZeta, 0,
            0, 0, 1;

    for (int i = 0; i < P.rows(); i++) {
        P.row(i) =  P.row(i) * Rx * Ry * Rz;
    }

    for (int i = 0; i < N.rows(); i++) {
        N.row(i) = N.row(i) * Rx * Ry * Rz;
    }

}

// Wendland weight function
double weightFunction(double r) {
    return pow(1-r/h, 4)*(4*r/h+1);
}

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

    bb_min *= 1.1;
    bb_max *= 1.1;

    // Bounding box dimensions
    Eigen::RowVector3d dim = bb_max - bb_min;

    // Grid spacing
    const double dx = dim[0] / (double)(resolutionX - 1);
    const double dy = dim[1] / (double)(resolutionY - 1);
    const double dz = dim[2] / (double)(resolutionZ - 1);

    // 3D positions of the grid points -- see slides or marching_cubes.h for ordering
    grid_points.resize(resolutionX * resolutionY * resolutionZ, 3);
    // Create each gridpoint
    for (unsigned int x = 0; x < resolutionX; ++x) {
        for (unsigned int y = 0; y < resolutionY; ++y) {
            for (unsigned int z = 0; z < resolutionZ; ++z) {
                // Linear index of the point at (x,y,z)
                int index = x + resolutionX * (y + resolutionY * z);
                // 3D point at (x,y,z)
                grid_points.row(index) = bb_min + Eigen::RowVector3d(x * dx, y * dy, z * dz);
            }
        }
    }

    grid_values.resize(resolutionX * resolutionY * resolutionZ);
    grid_values.setZero(resolutionX * resolutionY * resolutionZ);

}

// Function for explicitly evaluating the implicit function for a sphere of
// radius r centered at c : f(p) = ||p-c|| - r, where p = (x,y,z).
// This will NOT produce valid results for any mesh other than the given
// sphere.
// Replace this with your own function for evaluating the implicit function
// values at the grid points using MLS
void evaluateImplicitFunc() {

    // Bounding box dimensions
    Eigen::RowVector3d bb_min, bb_max;
    bb_min = P.colwise().minCoeff();
    bb_max = P.colwise().maxCoeff();
    Eigen::RowVector3d dim = bb_max - bb_min;

    double wendlandRadius = wendlandRadius_factor * sqrt(pow(dim[0], 2) + pow(dim[1], 2) + pow(dim[2], 2));
    
    // Polinom Degree

    int AColSize = 4;

    switch (polyDegree) {
        case 0:
            AColSize = 1;
        case 1:
            AColSize = 4;
        break;
        case 2:
            AColSize = 10;
        break;
    }

    // Evaluate sphere's signed distance function at each gridpoint.

    int counter = 0;

    for (unsigned int x = 0; x < resolutionX; ++x) {
        for (unsigned int y = 0; y < resolutionY; ++y) {
            for (unsigned int z = 0; z < resolutionZ; ++z) {

                // Linear index of the point at (x,y,z)

                int index = x + resolutionX * (y + resolutionY * z);

                // Current grid point

                Eigen::RowVector3d point1;
                point1(0) = grid_points(index, 0);
                point1(1) = grid_points(index, 1);
                point1(2) = grid_points(index, 2);

                // Initialization

                Eigen::MatrixXd A;
                Eigen::MatrixXd f;
                Eigen::MatrixXd w;

                A.resize(0, AColSize);
                f.resize(0, 1);
                w.resize(0, 0);

                // Searching closest neighbors

                int xRaw = floor(((abs(bb_min(0)) + grid_points(index, 0)) / dim[0])*resolutionXRaw);
                int yRaw = floor(((abs(bb_min(1)) + grid_points(index, 1)) / dim[1])*resolutionYRaw);
                int zRaw = floor(((abs(bb_min(2)) + grid_points(index, 2)) / dim[2])*resolutionZRaw);

                for (int x = xRaw - 1; x < xRaw + 2; x++) {

                    if (x < 0 || x >= resolutionXRaw) {
                        continue;
                    }

                    for (int y = yRaw - 1; y < yRaw + 2; y++) {

                        if (y < 0 || y >= resolutionYRaw) {
                            continue;
                        }

                        for (int z = zRaw - 1; z < zRaw + 2; z++) {

                            if (z < 0 || z >= resolutionZRaw) {
                                continue;
                            }

                            int indexRaw = x + resolutionXRaw * (y + resolutionYRaw * z);

                            for (int i : constrained_grid_points[indexRaw]) {

                                Eigen::RowVector3d point3;

                                point3(0) = constrained_points(i, 0);
                                point3(1) = constrained_points(i, 1);
                                point3(2) = constrained_points(i, 2);

                                Eigen::RowVector3d diff = point1 - point3;

                                double distance = sqrt(pow(diff(0), 2) + pow(diff(1), 2) + pow(diff(2), 2));

                                if (distance < wendlandRadius) {

                                    A.conservativeResize(A.rows() + 1, AColSize);
                                    f.conservativeResize(f.rows() + 1, 1);
                                    w.conservativeResize(w.rows() + 1, w.rows() + 1);

                                    // Update A

                                    A(A.rows() - 1, 0) = 1;
                                    if (polyDegree != 0) {
                                        for (int j = 0; j < 3; j++) {
                                            A(A.rows() - 1, j + 1) = constrained_points(i, j);
                                        }
                                        if (polyDegree == 2) {
                                            A(A.rows() - 1, 4) = constrained_points(i, 0)*constrained_points(i, 0);
                                            A(A.rows() - 1, 5) = constrained_points(i, 1)*constrained_points(i, 1);
                                            A(A.rows() - 1, 6) = constrained_points(i, 2)*constrained_points(i, 2);

                                            A(A.rows() - 1, 7) = constrained_points(i, 0)*constrained_points(i, 1);
                                            A(A.rows() - 1, 8) = constrained_points(i, 1)*constrained_points(i, 2);
                                            A(A.rows() - 1, 9) = constrained_points(i, 0)*constrained_points(i, 2);
                                        }
                                    }

                                    // Update f

                                    f(f.rows() - 1, 0) = constrained_values(i, 0);

                                    // Update w

                                    for (int j = 0; j < w.rows(); j++) {
                                        w(j, w.rows() - 1) = 0;
                                        w(w.rows() - 1, j) = 0;
                                    }
                                    w(w.rows() - 1, w.rows() - 1) = weightFunction(distance);

                                }

                            }
                        }
                    }
                }

                if (A.rows() != 0) {

                    // Solving the linear equations

                    Eigen::MatrixXd A_tag = w * A;
                    Eigen::MatrixXd f_tag = w * f;

                    Eigen::RowVectorXd c;
                    c.resize(1, AColSize);
                    c = A.colPivHouseholderQr().solve(f_tag).transpose();

                    // value of current grid point

                    Eigen::RowVectorXd point2;
                    point2.resize(1, AColSize);
                    point2(0) = 1;
                    if (polyDegree != 0) {
                        point2(1) = grid_points(index, 0);
                        point2(2) = grid_points(index, 1);
                        point2(3) = grid_points(index, 2);
                        if (polyDegree == 2) {
                            point2(4) = grid_points(index, 0)*grid_points(index, 0);
                            point2(5) = grid_points(index, 1)*grid_points(index, 1);
                            point2(6) = grid_points(index, 2)*grid_points(index, 2);

                            point2(7) = grid_points(index, 0)*grid_points(index, 1);
                            point2(8) = grid_points(index, 1)*grid_points(index, 2);
                            point2(9) = grid_points(index, 0)*grid_points(index, 2);
                        }
                    }

                    grid_values[index] = point2.dot(c);

                }
                else {

                    grid_values[index] = MAX_INT;

                }

                printf("%d of %d\n", ++counter, resolutionX*resolutionY*resolutionZ);

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

    for (unsigned int x = 0; x<resolutionZ; ++x) {
        for (unsigned int y = 0; y < resolutionY; ++y) {
            for (unsigned int z = 0; z < resolutionZ; ++z) {
                int index = x + resolutionX * (y + resolutionY * z);
                if (x < resolutionZ - 1) {
                    int index1 = (x + 1) + y * resolutionX + z * resolutionX * resolutionY;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (y < resolutionY - 1) {
                    int index1 = x + (y + 1) * resolutionX + z * resolutionX * resolutionY;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (z < resolutionZ - 1) {
                    int index1 = x + y * resolutionX + (z + 1) * resolutionX * resolutionY;
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

        // Bounding box dimensions
        Eigen::RowVector3d bb_min, bb_max;
        bb_min = P.colwise().minCoeff();
        bb_max = P.colwise().maxCoeff();
        Eigen::RowVector3d dim = bb_max - bb_min;

        double epsilon = epsilon_factor * sqrt(pow(dim[0], 2) + pow(dim[1], 2) + pow(dim[2], 2));

        constrained_points.resize(3 * P.rows(), 3);
        constrained_values.resize(3 * P.rows(), 1);

        Eigen::MatrixXd colors_per_face;
        colors_per_face.setZero(3 * P.rows(), 3);
        
        // Constrain points grid

        if (constrained_grid_points != NULL) {
            delete[] constrained_grid_points;
        }

        constrained_grid_points = new list<int>[resolutionXRaw * resolutionYRaw * resolutionZRaw];

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

            // add constrain point to raw grid

            int x = floor(((abs(bb_min(0)) + P(i, 0))/dim[0])*resolutionXRaw);
            if (x < 0) {
                x = 0;
            }
            else if (x >= resolutionXRaw) {
                x = resolutionXRaw - 1;
            }

            int y = floor(((abs(bb_min(1)) + P(i, 1))/dim[1])*resolutionYRaw);
            if (y < 0) {
                y = 0;
            }
            else if (y >= resolutionYRaw) {
                y = resolutionYRaw - 1;
            }

            int z = floor(((abs(bb_min(2)) + P(i, 2))/dim[2])*resolutionZRaw);
            if (z < 0) {
                z = 0;
            }
            else if (z >= resolutionZRaw) {
                z = resolutionZRaw - 1;
            }

            int index = x + resolutionXRaw * (y + resolutionYRaw * z);

            constrained_grid_points[index].push_back(3 * i);
            constrained_grid_points[index].push_back(3 * i + 1);
            constrained_grid_points[index].push_back(3 * i + 2);

        }

        // Add code for displaying all points, as above

        viewer.data().point_size = 11;
        viewer.data().add_points(constrained_points, colors_per_face);

    }

    if (key == '3') {
        // Show grid points with colored nodes and connected with lines
        viewer.data().clear();
        viewer.core.align_camera_center(P);

        // Make grid
        createGrid();

        // Evaluate implicit function
        evaluateImplicitFunc();

        // get grid lines
        getLines();

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

       viewer.data().add_edges(grid_lines.block(0, 0, grid_lines.rows(), 3), grid_lines.block(0, 3, grid_lines.rows(), 3), Eigen::RowVector3d(0.8, 0.8, 0.8));

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
        igl::copyleft::marching_cubes(grid_values, grid_points, resolutionX, resolutionY, resolutionZ, V, F);
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
        ImGui::InputInt("Resolution-X", &resolutionX, 0, 0);
        ImGui::InputInt("Resolution-Y", &resolutionY, 0, 0);
        ImGui::InputInt("Resolution-Z", &resolutionZ, 0, 0);

        ImGui::InputDouble("Wendland Radius", &wendlandRadius_factor, 0, 0);
        ImGui::InputDouble("Epsilon", &epsilon_factor, 0, 0);
        ImGui::InputInt("Polinom Degree", &polyDegree, 0, 0);
        ImGui::InputDouble("Smoothing Factor", &h, 0, 0);

        ImGui::InputInt("Object", &object, 0, 0);
        if (ImGui::Button("Load Object", ImVec2(-1, 0)))
        {
            switch (object) {
            case 0:
                igl::readOFF("../data/Test.off", P, F, N);
                break;
            case 1:
                igl::readOFF("../data/sphere.off", P, F, N);
                break;
            case 2:
                igl::readOFF("../data/cat.off", P, F, N);
                break;
            case 3:
                igl::readOFF("../data/luigi.off", P, F, N);
                break;
            case 4:
                igl::readOFF("../data/hound.off", P, F, N);
                break;
            case 5:
                igl::readOFF("../data/horse.off", P, F, N);
                break;
            case 6:
                igl::readOFF("../data/bunny-500.off", P, F, N);
                break;
            case 7:
                igl::readOFF("../data/bunny-1000.off", P, F, N);
                break;
            }

            rotate();

            callback_key_down(viewer, '1', 0);
        }

        ImGui::InputInt("Resolution-X-index", &resolutionXRaw, 0, 0);
        ImGui::InputInt("Resolution-Y-index", &resolutionYRaw, 0, 0);
        ImGui::InputInt("Resolution-Z-index", &resolutionZRaw, 0, 0);

        ImGui::InputDouble("Theta", &theta, 0, 0);
        ImGui::InputDouble("Phi", &phi, 0, 0);
        ImGui::InputDouble("Zeta", &zeta, 0, 0);

      }

    };

    viewer.launch();
}
