/** @file lrmb.cpp
 *
 * Given two triangle meshes, perform the boolean operation and return
 * a third mesh.
 */

#include "lrmb.h"
#include "booleans.h"

extern "C" int
mesh_boolean(
	int **faces, int *num_faces, point_t **vertices, int *num_vertices,
	const int *f1, int nf1, const point_t *v1, int nv1,
	int op,
	const int *f2, int nf2, const point_t *v2, int nv2
	)
{
    if (!faces || !num_faces || !vertices || !num_vertices) return -1;
    if (!f1|| !nf1 || !v1 || !nv1) return -1;
    if (!f2|| !nf2 || !v2 || !nv2) return -1;

    // Prepare the IRMB containers
    std::vector<double> in_coords;
    std::vector<uint> in_tris;
    std::vector<uint> in_labels;

    // Load first mesh data
    for (int i = 0; i < nv1; i++) {
	in_coords.push_back(v1[i][0]);
	in_coords.push_back(v1[i][1]);
	in_coords.push_back(v1[i][2]);
    }
    for (int i = 0; i < nf1*3; i++) {
	in_tris.push_back(f1[i]);
    }
    for (int i = 0; i < nf1; i++)
	in_labels.push_back(0);

    // Second mesh's vertex indices will be
    // offset in the in_coord array
    size_t vert_offset = (size_t)(in_coords.size() / 3);

    // Load second mesh data
    for (int i = 0; i < nv2; i++) {
	in_coords.push_back(v2[i][0]);
	in_coords.push_back(v2[i][1]);
	in_coords.push_back(v2[i][2]);
    }
    for (int i = 0; i < nf2*3; i++) {
	in_tris.push_back(f2[i] + vert_offset);
    }
    for (int i = 0; i < nf2; i++)
	in_labels.push_back(1);

    // Translate boolean
    BoolOp o_op;
    switch (op) {
	case 1:
	    o_op = SUBTRACTION;
	    break;
	case 2:
	    o_op = INTERSECTION;
	    break;
	default:
	    o_op = UNION;
    }

    // Run the pipeline
    std::vector<double> bool_coords;
    std::vector<uint> bool_tris;
    std::vector<std::bitset<32>> bool_labels;
    booleanPipeline(in_coords, in_tris, in_labels, o_op, bool_coords, bool_tris, bool_labels);


    return 0;
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8

