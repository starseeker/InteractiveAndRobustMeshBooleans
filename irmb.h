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

#ifndef IRMB_DEFINES_H
#define IRMB_DEFINES_H

#ifndef IRMB_EXPORT
#  if defined(IRMB_DLL_EXPORTS) && defined(IRMB_DLL_IMPORTS)
#    error "Only IRMB_DLL_EXPORTS or IRMB_DLL_IMPORTS can be defined, not both."
#  elif defined(IRMB_DLL_EXPORTS)
#    define IRMB_EXPORT COMPILER_DLLEXPORT
#  elif defined(IRMB_DLL_IMPORTS)
#    define IRMB_EXPORT COMPILER_DLLIMPORT
#  else
#    define IRMB_EXPORT
#  endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

IRMB_EXPORT extern long
bool_meshes(
    double **o_coords, int *o_clen, unsigned int **o_tris, int *o_tricnt,
    int b_op,
    double *a_coords, int a_clen, unsigned int *a_tris, int a_tricnt,
    double *b_coords, int b_clen, unsigned int *b_tris, int b_tricnt
    );

#ifdef __cplusplus
}
#endif

#endif /* IRMB_DEFINES_H */

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
