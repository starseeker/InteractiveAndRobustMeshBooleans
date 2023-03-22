/*****************************************************************************************
 *              MIT License                                                              *
 *                                                                                       *
 * Copyright (c) 2022 G. Cherchi, F. Pellacini, M. Attene and M. Livesu                  *
 *                                                                                       *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this  *
 * software and associated documentation files (the "Software"), to deal in the Software *
 * without restriction, including without limitation the rights to use, copy, modify,    *
 * merge, publish, distribute, sublicense, and/or sell copies of the Software, and to    *
 * permit persons to whom the Software is furnished to do so, subject to the following   *
 * conditions:                                                                           *
 *                                                                                       *
 * The above copyright notice and this permission notice shall be included in all copies *
 * or substantial portions of the Software.                                              *
 *                                                                                       *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,   *
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A         *
 * PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT    *
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION     *
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE        *
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                                *
 *                                                                                       *
 * Authors:                                                                              *
 *      Gianmarco Cherchi (g.cherchi@unica.it)                                           *
 *      https://www.gianmarcocherchi.com                                                 *
 *                                                                                       *
 *      Fabio Pellacini (fabio.pellacini@uniroma1.it)                                    *
 *      https://pellacini.di.uniroma1.it                                                 *
 *                                                                                       *
 *      Marco Attene (marco.attene@ge.imati.cnr.it)                                      *
 *      https://www.cnr.it/en/people/marco.attene/                                       *
 *                                                                                       *
 *      Marco Livesu (marco.livesu@ge.imati.cnr.it)                                      *
 *      http://pers.ge.imati.cnr.it/livesu/                                              *
 *                                                                                       *
 * ***************************************************************************************/

#include <cinolib/meshes/trimesh.h>
#include <cinolib/find_intersections.h>
#include "booleans.h"
#include "irmb.h"

static bool
negative_mesh(double *coords, int clen, unsigned int *tris, int tricnt)
{
    std::vector<double> cvec;
    for (int i = 0; i < clen/3; i++) {
	cvec.push_back(coords[3*i]);
	cvec.push_back(coords[3*i+1]);
	cvec.push_back(coords[3*i+2]);
	std::cout << "V(" << i << "):" << coords[3*i] << "," << coords[3*i+1] << "," << coords[3*i+2] << "\n";
    }
    std::vector<unsigned int> tvec;
    for (int i = 0; i < tricnt; i++) {
	tvec.push_back(tris[3*i]);
	tvec.push_back(tris[3*i+1]);
	tvec.push_back(tris[3*i+2]);
	std::cout << "T(" << i << "):" << tris[3*i] << "," << tris[3*i+1] << "," << tris[3*i+2] << "\n";
    }

    cinolib::Trimesh<> m(cvec, tvec);

    // Global orientaiton
    // WARNING: this is correct only if the mesh has a single connected component
    cinolib::vec3d  O   = m.bbox().min - cinolib::vec3d(1,0,0); // surely external
    double vol = 0;
    for(uint pid=0; pid<m.num_polys(); ++pid)
    {
	// sum signed volumes
	vol += tet_volume(m.poly_vert(pid,0),
		m.poly_vert(pid,1),
		m.poly_vert(pid,2), O);
    }
    // if faces are CCW volume is negative
    if (vol<0) {
	std::cout << "- CCW volume - negative mesh!" << std::endl;
	return true;
    }

    return false;
}

static bool
mesh_valid(double *coords, int clen, unsigned int *tris, int tricnt)
{
    if (!coords || !clen || !tris || !tricnt)
	return false;

    std::vector<double> cvec;
    for (int i = 0; i < clen/3; i++) {
	cvec.push_back(coords[3*i]);
	cvec.push_back(coords[3*i+1]);
	cvec.push_back(coords[3*i+2]);
	std::cout << "V(" << i << "):" << coords[3*i] << "," << coords[3*i+1] << "," << coords[3*i+2] << "\n";
    }
    std::vector<unsigned int> tvec;
    for (int i = 0; i < tricnt; i++) {
	tvec.push_back(tris[3*i]);
	tvec.push_back(tris[3*i+1]);
	tvec.push_back(tris[3*i+2]);
	std::cout << "T(" << i << "):" << tris[3*i] << "," << tris[3*i+1] << "," << tris[3*i+2] << "\n";
    }

    cinolib::Trimesh<> m(cvec, tvec);

    // CHECK if mesh is MANIFOLD
    for(uint vid=0; vid<m.num_verts(); ++vid) {
	if(!m.vert_is_manifold(vid)) {
	    std::cout << "- your input is NOT manifold!" << std::endl;
	    std::cout << "- CHECK FAILED!" << std::endl;
	    return false;
	}
    }

    // CHECK triangle ORIENTATION
    {
	// Local orientation
	std::vector<uint> visited(m.num_polys(),false);
	for(uint pid=0; pid<m.num_polys(); ++pid)
	{
	    if(visited[pid]) continue;
	    std::queue<uint> q;
	    q.push(pid);
	    visited[pid] = true;
	    while(!q.empty())
	    {
		uint pid = q.front();
		q.pop();
		for(uint nbr : m.adj_p2p(pid))
		{
		    uint eid = m.edge_shared(pid, nbr);
		    if(m.edge_is_CCW(eid,pid)==m.edge_is_CCW(eid,nbr)) {
			std::cout << "- your mesh local orientation is NOT correct!" << std::endl;
			std::cout << "- CHECK FAILED!" << std::endl;
			return false;
		    }
		    if(!visited[nbr])
		    {
			visited[nbr] = true;
			q.push(nbr);
		    }
		}
	    }
	}
    }

    // CHECK self INTERSECTIONS
    std::set<cinolib::ipair> intersections;
    cinolib::find_intersections(m, intersections);

    if(!intersections.empty())
    {
	std::cout << "- your input is NOT self-intersections free!" << std::endl;
	std::cout << "- CHECK FAILED!" << std::endl;
	return false;
    }

    return true;
}

extern "C" long
bool_meshes(
	double **o_coords, int *o_clen, unsigned int **o_tris, int *o_tricnt,
	int b_op,
	double *a_coords, int a_clen, unsigned int *a_tris, int a_tricnt,
	double *b_coords, int b_clen, unsigned int *b_tris, int b_tricnt
	)
{
    // If the mesh is negative, flip the triangles
    if (negative_mesh(a_coords, a_clen, a_tris, a_tricnt)) {
	for (int i = 0; i < a_tricnt; i++) {
	    int tmp_tri = a_tris[3*i+1];
	    a_tris[3*i+1] = a_tris[3*i+2];
	    a_tris[3*i+2] = tmp_tri;
	}
    }
    if (negative_mesh(a_coords, a_clen, a_tris, a_tricnt)) {
	std::cout << "a flipping failed!\n";
    }
    if (negative_mesh(b_coords, b_clen, b_tris, b_tricnt)) {
	for (int i = 0; i < b_tricnt; i++) {
	    int tmp_tri = b_tris[3*i+1];
	    b_tris[3*i+1] = b_tris[3*i+2];
	    b_tris[3*i+2] = tmp_tri;
	}
    }
    if (negative_mesh(b_coords, b_clen, b_tris, b_tricnt)) {
	std::cout << "b flipping failed!\n";
    }

    // First step, check the input meshes.  If they don't satisfy the criteria, don't proceed
    if (!mesh_valid(a_coords, a_clen, a_tris, a_tricnt))
	return -1;
    if (!mesh_valid(b_coords, b_clen, b_tris, b_tricnt))
	return -1;

    // Decode boolean operator
    BoolOp op;
    switch (b_op) {
	case 1:
	    op = SUBTRACTION;
	    break;
	case 2:
	    op = INTERSECTION;
	    break;
	case 3:
	    op = XOR;
	    break;
	default:
	    op = UNION;
	    break;
    }

    // Pack the data for boolean pipeline processing
    std::vector<double> in_coords;
    std::vector<uint> in_tris;
    std::vector<uint> in_labels;

    for (int i = 0; i < a_clen; i++) in_coords.push_back(a_coords[i]);
    size_t offset = static_cast<size_t>(in_coords.size() / 3);
    for (int i = 0; i < a_tricnt*3; i++) in_tris.push_back(a_tris[i]);
    for (int i = 0; i < a_tricnt; i++) in_labels.push_back(0);

    for (int i = 0; i < b_clen; i++) in_coords.push_back(b_coords[i]);
    for (int i = 0; i < b_tricnt*3; i++) in_tris.push_back(b_tris[i]+offset);
    for (int i = 0; i < b_tricnt; i++) in_labels.push_back(1);

    // Perform the boolean pipeline
    std::vector<double> bool_coords;
    std::vector<uint> bool_tris;
    std::vector<std::bitset<NBIT>> bool_labels;
    booleanPipeline(in_coords, in_tris, in_labels, op, bool_coords, bool_tris, bool_labels);

    // If we don't have the output targets, just report the final triangle count
    if (!o_coords || !o_clen || !o_tris || !o_tricnt)
	return static_cast<long>(bool_tris.size() / 3);

    // If we don't have valid output, we have nothing
    if (!bool_coords.size() || !bool_tris.size())
	return 0;

    // Success - populate the output containers
    (*o_coords) = (double *)calloc(bool_coords.size(), sizeof(double));
    (*o_tris) = (unsigned int *)calloc(bool_tris.size(), sizeof(unsigned int));
    for (size_t i = 0; i < bool_coords.size(); i++) (*o_coords)[i] = bool_coords[i];
    for (size_t i = 0; i < bool_tris.size(); i++) (*o_tris)[i] = bool_tris[i];
    *o_clen = static_cast<int>(bool_coords.size()/3);
    *o_tricnt = static_cast<int>(bool_tris.size()/3);

    return static_cast<long>(bool_tris.size()/3);
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8

