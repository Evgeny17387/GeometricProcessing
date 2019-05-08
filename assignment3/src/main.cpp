#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <igl/local_basis.h>
#include <igl/grad.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>


/*** insert any necessary libigl headers here ***/
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/lscm.h>
#include <igl/adjacency_matrix.h>
#include <igl/sum.h>
#include <igl/diag.h>
#include <igl/speye.h>
#include <igl/repdiag.h>
#include <igl/cat.h>
#include <igl/dijkstra.h>
#include <igl/doublearea.h>

#include <cmath>

using namespace std;
using namespace Eigen;
using Viewer = igl::opengl::glfw::Viewer;

Viewer viewer;

// vertex array, #V x3
Eigen::MatrixXd V;

// face array, #F x3
Eigen::MatrixXi F;

// UV coordinates, #V x2
Eigen::MatrixXd UV;

bool showingUV = false;
bool freeBoundary = false;
double TextureResolution = 10;
igl::opengl::ViewerCore temp3D;
igl::opengl::ViewerCore temp2D;

bool isUV = false;
bool isDeformation = false;
int deformationMode = 0;

void Redraw()
{
	viewer.data().clear();

	if (!showingUV)
	{
		viewer.data().set_mesh(V, F);
		viewer.data().set_face_based(false);

        if(UV.size() != 0)
        {
          viewer.data().set_uv(TextureResolution*UV);
          viewer.data().show_texture = true;
        }
	}
    else
	{
		viewer.data().show_texture = false;
		viewer.data().set_mesh(UV, F);
	}

    if (isUV && isDeformation) {

        Eigen::MatrixXd colors_per_face;

        colors_per_face.resize(F.rows(), 3);
        colors_per_face.setZero();

        double maxDeformation = 0;

        for (int i = 0; i < F.rows(); i++) {

            double face_color;

            double deformation;

            if (deformationMode == 0) {

                double edge3D1 = (V.row(F(i, 0)) - V.row(F(i, 1))).norm();
                double edge3D2 = (V.row(F(i, 1)) - V.row(F(i, 2))).norm();
                double edge3D3 = (V.row(F(i, 2)) - V.row(F(i, 0))).norm();

                double edge3DSum = edge3D1 + edge3D2 + edge3D3;

                double edge2D1 = (UV.row(F(i, 0)) - UV.row(F(i, 1))).norm();
                double edge2D2 = (UV.row(F(i, 1)) - UV.row(F(i, 2))).norm();
                double edge2D3 = (UV.row(F(i, 2)) - UV.row(F(i, 0))).norm();

                double edge2DSum = edge2D1 + edge2D2 + edge2D3;

                deformation = abs(edge3DSum - edge2DSum);

            }
            else if (deformationMode == 1)
            {

                deformation = 0;

            }

            if (deformation > maxDeformation) {
                maxDeformation = deformation;
            }

            colors_per_face(i, 0) = deformation;

        }

        for (int i = 0; i < F.rows(); i++) {

            double deformationFactor = colors_per_face(i, 0) / maxDeformation;

            if (deformationFactor > 0.75) {

                colors_per_face(i, 0) = 0.75 + 1 - deformationFactor;
                colors_per_face(i, 1) = 0;
                colors_per_face(i, 2) = 0;

            }
            else
            {

                colors_per_face(i, 0) = 1;
                colors_per_face(i, 1) = 1 - deformationFactor;
                colors_per_face(i, 2) = 1 - deformationFactor;

            }

        }

        viewer.data().set_colors(colors_per_face);

    }

}

bool callback_mouse_move(Viewer &viewer, int mouse_x, int mouse_y)
{
	if (showingUV)
		viewer.mouse_mode = igl::opengl::glfw::Viewer::MouseMode::Translation;
	return false;
}

static void computeSurfaceGradientMatrix(SparseMatrix<double> & D1, SparseMatrix<double> & D2)
{
	MatrixXd F1, F2, F3;
	SparseMatrix<double> DD, Dx, Dy, Dz;

	igl::local_basis(V, F, F1, F2, F3);
	igl::grad(V, F, DD);

	Dx = DD.topLeftCorner(F.rows(), V.rows());
	Dy = DD.block(F.rows(), 0, F.rows(), V.rows());
	Dz = DD.bottomRightCorner(F.rows(), V.rows());

	D1 = F1.col(0).asDiagonal()*Dx + F1.col(1).asDiagonal()*Dy + F1.col(2).asDiagonal()*Dz;
	D2 = F2.col(0).asDiagonal()*Dx + F2.col(1).asDiagonal()*Dy + F2.col(2).asDiagonal()*Dz;
}
static inline void SSVD2x2(const Eigen::Matrix2d& J, Eigen::Matrix2d& U, Eigen::Matrix2d& S, Eigen::Matrix2d& V)
{
	double e = (J(0) + J(3))*0.5;
	double f = (J(0) - J(3))*0.5;
	double g = (J(1) + J(2))*0.5;
	double h = (J(1) - J(2))*0.5;
	double q = sqrt((e*e) + (h*h));
	double r = sqrt((f*f) + (g*g));
	double a1 = atan2(g, f);
	double a2 = atan2(h, e);
	double rho = (a2 - a1)*0.5;
	double phi = (a2 + a1)*0.5;

	S(0) = q + r;
	S(1) = 0;
	S(2) = 0;
	S(3) = q - r;

	double c = cos(phi);
	double s = sin(phi);
	U(0) = c;
	U(1) = s;
	U(2) = -s;
	U(3) = c;

	c = cos(rho);
	s = sin(rho);
	V(0) = c;
	V(1) = -s;
	V(2) = s;
	V(3) = c;
}

void ConvertConstraintsToMatrixForm(VectorXi indices, MatrixXd positions, SparseMatrix<double> &C, VectorXd &d)
{
	// Convert the list of fixed indices and their fixed positions to a linear system
	// Hint: The matrix C should contain only one non-zero element per row and d should contain the positions in the correct order.

    C.resize(2 * indices.rows(), 2 * V.rows());
    d.resize(2 * indices.rows(), 1);

    C.setZero();

    for (int i = 0; i < indices.rows(); i++) {

        C.coeffRef(i, indices(i)) = 1;
        C.coeffRef(indices.rows() + i, V.rows() + indices(i)) = 1;

    }

    for (int i = 0; i < indices.rows(); i++) {
        d(i, 0) = positions(i, 0);
    }

    for (int i = 0; i < indices.rows(); i++) {
        d(indices.rows() + i, 0) = positions(i, 1);
    }

}

void printVectorOfVectors(std::vector<std::vector<int>> VF) {
    int i = -1;
    for (auto vector1 : VF) {
        printf("%d:", ++i);
        for (auto vector2 : vector1) {
            printf(" %d", vector2);
        }
        printf("\n");
    }
    printf("\n");
}

void printVectorXi(VectorXi vectorXi) {
    for (int i = 0; i < vectorXi.rows(); i++) {
        printf("%d\n", vectorXi(i));
    }
}

void printVectorXd(VectorXd MatrixXd) {
    for (int i = 0; i < MatrixXd.rows(); i++) {
        printf("%f\n", MatrixXd(i));
    }
}

void printMatrixXd(MatrixXd matrixXd) {
    for (int i = 0; i < matrixXd.rows(); i++) {
        for (int j = 0; j < matrixXd.cols(); j++) {
            printf("%f ", matrixXd(i, j));
        }
        printf("\n");
    }
}

void computeParameterization(int type)
{
	VectorXi fixed_UV_indices;
	MatrixXd fixed_UV_positions;

	SparseMatrix<double> A;
	VectorXd b;
	SparseMatrix<double> C;
	VectorXd d;

    // Find the indices of the boundary vertices of the mesh and put them in fixed_UV_indices
	if (!freeBoundary)
	{
		// The boundary vertices should be fixed to positions on the unit disc. Find these position and
		// save them in the #V x 2 matrix fixed_UV_position.

        igl::boundary_loop(F, fixed_UV_indices);

        igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);

    }
	else
	{
		// Fix two UV vertices. This should be done in an intelligent way. Hint: The two fixed vertices should be the two most distant one on the mesh.

        int V1;
        int V2;

        int source = 0;

        set<int> targets;

        vector<vector<int>> VV;
        igl::adjacency_list(F, VV);

        VectorXd min_distance;
        VectorXd previous;

        igl::dijkstra(source, targets, VV, min_distance, previous);

        int max_distance = min_distance(0);
        source = 0;
        for (int i = 1; i < min_distance.size(); i++) {
            if (max_distance < min_distance(i)) {
                max_distance = min_distance(i);
                source = i;
            }
        }

        V1 = source;

        igl::dijkstra(source, targets, VV, min_distance, previous);

        max_distance = min_distance(0);
        source = 0;
        for (int i = 1; i < min_distance.size(); i++) {
            if (max_distance < min_distance(i)) {
                max_distance = min_distance(i);
                source = i;
            }
        }

        V2 = source;

        fixed_UV_indices.resize(2, 1);
        fixed_UV_indices(0, 0) = V1;
        fixed_UV_indices(1, 0) = V2;

        fixed_UV_positions.resize(2, 3);
        fixed_UV_positions.row(0) = V.row(V1);
        fixed_UV_positions.row(1) = V.row(V2);

    }

    if ((type == '1') || (type == '2') || (type == '3') || (type == '4')) {

        ConvertConstraintsToMatrixForm(fixed_UV_indices, fixed_UV_positions, C, d);

//        printVectorXi(fixed_UV_indices);

//    printMatrixXd(C);
//    printf("\n");

//    printVectorXd(d);
//    printf("\n");

    }

    // Find the linear system for the parameterization (1- Tutte, 2- Harmonic, 3- LSCM, 4- ARAP)
	// and put it in the matrix A.
	// The dimensions of A should be 2#V x 2#V.
	if (type == '1') {
        
        // Add your code for computing uniform Laplacian for Tutte parameterization
        // Hint: use the adjacency matrix of the mesh

        SparseMatrix<double> B;
        igl::adjacency_matrix(F, B);

        SparseVector<double> Bsum;
        igl::sum(B, 1, Bsum);

        SparseMatrix<double> D;
        igl::diag(Bsum, D);

        SparseMatrix<double> DInverse;
        DInverse.resize(D.rows(), D.cols());
        DInverse.setZero();

        for (int i = 0; i < D.rows(); i++) {
            DInverse.coeffRef(i, i) = 1 / D.coeffRef(i, i);
        }

//        printMatrixXd(DInverse);

        SparseMatrix<double> Identity;
        Identity.resize(D.rows(), D.cols());
        Identity.setIdentity();

//        printMatrixXd(Identity);

        SparseMatrix<double> L;
        L = Identity - DInverse * B;

//        printMatrixXd(L);

        SparseMatrix<double> Up;
        SparseMatrix<double> Down;

        SparseMatrix<double> Temp;
        Temp.resize(L.rows(), L.cols());
        Temp.setZero();

        igl::cat(2, L, Temp, Up);
        igl::cat(2, Temp, L, Down);
        igl::cat(1, Up, Down, A);

//        printMatrixXd(A);

        b.resize(2 * V.rows(), 1);
        b.setZero();

    }

	if (type == '2') {
		// Add your code for computing cotangent Laplacian for Harmonic parameterization
		// Use can use a function "cotmatrix" from libIGL, but ~~~~***READ THE DOCUMENTATION***~~~~

        SparseMatrix<double> L;

        igl::cotmatrix(V, F, L);

        SparseMatrix<double> Up;
        SparseMatrix<double> Down;

        SparseMatrix<double> Temp;
        Temp.resize(L.rows(), L.cols());
        Temp.setZero();

        SparseMatrix<double> LMinus = -L;

        igl::cat(2, L, Temp, Up);
        igl::cat(2, Temp, L, Down);
        igl::cat(1, Up, Down, A);

        b.resize(2 * V.rows(), 1);
        b.setZero();

    }

	if (type == '3') {
		// Add your code for computing the system for LSCM parameterization
		// Note that the libIGL implementation is different than what taught in the tutorial! Do not rely on it!!

        // Surface gradients

        SparseMatrix<double> D1;
        SparseMatrix<double> D2;

        computeSurfaceGradientMatrix(D1, D2);

//        printf("D1, D2\n");
//        printf("%d, %d\n", D1.rows(), D1.cols());
//        printf("%d, %d\n", D2.rows(), D2.cols());
//        printf("\n");

        // Triangles area

        MatrixXd AreaColumn;
        SparseMatrix<double> AreaMatrix;

        igl::doublearea(V, F, AreaColumn);

//        printf("AreaColumn\n");
//        printf("%d, %d\n", AreaColumn.rows(), AreaColumn.cols());
//        printf("\n");

        AreaMatrix.resize(AreaColumn.rows(), AreaColumn.rows());
        AreaMatrix.setZero();

        for (int i = 0; i < AreaMatrix.rows(); i++) {
            AreaMatrix.coeffRef(i, i) = AreaColumn(i, 0);
        }

//        printf("AreaMatix\n");
//        printMatrixXd(AreaMatrix);
//        printf("\n");

        // A

        SparseMatrix<double> a;
        SparseMatrix<double> Temp;

        SparseMatrix<double> Up;
        SparseMatrix<double> Down;

        a = D1.transpose()*AreaMatrix*D1 + D2.transpose()*AreaMatrix*D2;
        
        Temp.resize(a.rows(), a.cols());
        Temp.setZero();
            
        igl::cat(2, a, Temp, Up);
        igl::cat(2, Temp, a, Down);
        igl::cat(1, Up, Down, A);

//        printf("A\n");
//        printMatrixXd(A);
//        printf("\n");

        // b

        b.resize(2 * V.rows(), 1);
        b.setZero();

    }

	if (type == '4') {
		// Add your code for computing ARAP system and right-hand side
		// Implement a function that computes the local step first
		// Then construct the matrix with the given rotation matrices
	}

    if ((type == '1') || (type == '2') || (type == '3')) {

        // Solve the linear system.
        // Construct the system as discussed in class and the assignment sheet
        // Use igl::cat to concatenate matrices

        SparseMatrix<double> Up;
        SparseMatrix<double> Down;
        SparseMatrix<double> Final;

        SparseMatrix<double> CTranspose = C.transpose();

        SparseMatrix<double> Temp;
        Temp.resize(C.rows(), C.rows());
        Temp.setZero();

        igl::cat(2, A, CTranspose, Up);
        igl::cat(2, C, Temp, Down);
        igl::cat(1, Up, Down, Final);

        VectorXd bd;

        igl::cat(1, b, d, bd);

        // Use Eigen::SparseLU to solve the system. Refer to tutorial 3 for more detail

//        printMatrixXd(Final);
//        printf("\n");

        VectorXd x(bd.rows());
        SparseLU<SparseMatrix<double>> solver;
        solver.analyzePattern(Final);
        solver.factorize(Final);
        x = solver.solve(bd);

//        printVectorXd(x);
//        printf("\n");

        // The solver will output a vector
        UV.resize(V.rows(), 2);

        for (int i = 0; i < V.rows(); i++) {
            UV(i, 0) = x(i);
            UV(i, 1) = x(V.rows() + i);
        }

        isUV = true;

    }

}

bool callback_key_pressed(Viewer &viewer, unsigned char key, int modifiers) {
	switch (key) {
    case '0':
    case '1':
	case '2':
	case '3':
	case '4':
		computeParameterization(key);
		break;
	case '5':
			// Add your code for detecting and displaying flipped triangles in the
			// UV domain here
		break;
	case '+':
		TextureResolution /= 2;
		break;
	case '-':
		TextureResolution *= 2;
		break;
	case ' ': // space bar -  switches view between mesh and parameterization
    if(showingUV)
    {
      temp2D = viewer.core;
      viewer.core = temp3D;
      showingUV = false;
    }
    else
    {
      if(UV.rows() > 0)
      {
        temp3D = viewer.core;
        viewer.core = temp2D;
        showingUV = true;
      }
      else { std::cout << "ERROR ! No valid parameterization\n"; }
    }
    break;
	}
	Redraw();
	return true;
}

bool load_mesh(string filename)
{
	igl::read_triangle_mesh(filename,V,F);
	Redraw();
	viewer.core.align_camera_center(V);
	showingUV = false;

	return true;
}

bool callback_init(Viewer &viewer)
{
	temp3D = viewer.core;
	temp2D = viewer.core;
	temp2D.orthographic = true;

	return false;
}

int main(int argc,char *argv[]) {
	if(argc != 2) {
		cout << "Usage ex3_bin <mesh.off/obj>" << endl;
//        load_mesh("../data/Test.off");
//		    load_mesh("../data/cathead.obj");
//        load_mesh("../data/hemisphere.off");
        load_mesh("../data/Octo_cut2.obj");
//        load_mesh("../data/hemisphere_non_convex_boundary.obj");

    }
	else
	{
		// Read points and normals
		load_mesh(argv[1]);
	}

	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);

	menu.callback_draw_viewer_menu = [&]()
	{
		// Draw parent menu content
		menu.draw_viewer_menu();

		// Add new group
		if (ImGui::CollapsingHeader("Parmaterization", ImGuiTreeNodeFlags_DefaultOpen))
		{
			// Expose variable directly ...
			ImGui::Checkbox("Free boundary", &freeBoundary);

            ImGui::Checkbox("Deformation", &isDeformation);

            ImGui::InputInt("Deformation Mode", &deformationMode);

			// TODO: Add more parameters to tweak here...
		}
	};

	viewer.callback_key_pressed = callback_key_pressed;
	viewer.callback_mouse_move = callback_mouse_move;
	viewer.callback_init = callback_init;

	viewer.launch();
}
