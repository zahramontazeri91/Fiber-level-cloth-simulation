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
#include <random>

class Rng {
	std::mt19937_64 random_engine;
	
public:
	Rng(std::random_device &rd): random_engine(rd()) {}

	Rng(int seed): random_engine(seed) {}

	int index;

	float rand(float r_min, float r_max) {
		std::uniform_real_distribution<float> random_distribution(r_min, r_max);
		return random_distribution(random_engine);
	}
};