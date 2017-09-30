/*
Copyright (c) 2015 to 2016 by Cornell University and The Regents Of
The University Of California. All Rights Reserved.

Permission to use this Procedural Yarn Fitting and Generation Tool (the "Work")
and its associated copyrights solely for educational, research and non-profit
purposes, without fee is hereby granted, provided that the user agrees as
follows:

Those desiring to incorporate the Work into commercial products or use Work and
its associated copyrights for commercial purposes should contact the Center for
Technology Licensing at Cornell University at

395 Pine Tree Road, Suite 310, Ithaca, NY 14850;
email: ctl-connect@cornell.edu;
Tel: 607-254-4698;
FAX: 607-254-5454

for a commercial license.

IN NO EVENT SHALL CORNELL UNIVERSITY ("CORNELL") OR THE UNIVERSITY OF
CALIFORNIA ("UC") BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL,
INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF
THE USE OF THE WORK AND ITS ASSOCIATED COPYRIGHTS, EVEN IF CORNELL OR UC MAY
HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

THE WORK PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND NEITHER CORNELL NOR UC HAS
ANY OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
MODIFICATIONS. CORNELL AND UC MAKE NO REPRESENTATIONS AND EXTEND NO WARRANTIES
OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR
THAT THE USE OF WORK AND ITS ASSOCIATED COPYRIGHTS WILL NOT INFRINGE ANY PATENT,
TRADEMARK OR OTHER RIGHTS.
*/

#pragma once 
#include "Util.h"

struct CurvePoint {
	vec3f curve_point_vertex;
	vec3f curve_point_normal;
	vec2f curve_point_uv;
};

struct Curve {
	/* WEFT == U, WARP == V */
	typedef enum {WEFT, WARP} woven_t; 
	std::vector<bool>  curve_on_top;
	std::vector<vec3f> curve_vertices;	
	std::vector<vec3f> curve_normals;
	std::vector<vec2f> curve_uvs;
	woven_t			   curve_type;	
	std::vector<float> curve_length_accumulated;
	float			   curve_length;
};

typedef std::vector<Curve> Curves;