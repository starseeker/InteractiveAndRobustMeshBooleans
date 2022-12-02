extern "C" {
#include "vmath.h"
}

extern "C" int
mesh_boolean(
	int **faces, int *num_faces, point_t **vertices, int *num_vertices,
	const int *f1, int nf1, const point_t *v1, int nv1,
	int op,
	const int *f2, int nf2, const point_t *v2, int nv2
	);

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8

