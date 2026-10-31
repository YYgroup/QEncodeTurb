/*
 Authors: Nikita Gourianov and Dr. Michael Lubasch
 Date: $LastChangedDate: 2018 $
 (c) University of Oxford 2018
 */

 /*! \page softwarelic Software Licence
 *
 * TENSOR THEORY NETWORK LIBRARY © 2014, The University of Oxford (the "Software")
 *
 * The Software remains the property of the University of Oxford ("the University").
 *
 * The Software is distributed "AS IS" under this Licence solely for non-commercial use in the hope that it will be useful,
 * but in order that the University as a charitable foundation protects its assets for the benefit of its educational and
 * research purposes, the University makes clear that no condition is made or to be implied, nor is any warranty given or
 * to be implied, as to the accuracy of the Software, or that it will be suitable for any particular purpose or for use
 * under any specific conditions. Furthermore, the University disclaims all responsibility for the use which is made of
 * the Software. It further disclaims any liability for the outcomes arising from using the Software.
 *
 * The Licensee agrees to indemnify the University and hold the University harmless from and against any and all claims,
 * damages and liabilities asserted by third parties (including claims for negligence), which arise directly or indirectly
 * from the use of the Software or the sale of any products based on the Software.
 *
 * No part of the Software may be reproduced, modified, transmitted or transferred in any form or by any means, electronic
 * or mechanical, without the express permission of the University. The permission of the University is not required if the
 * said reproduction, modification, transmission or transference is done without financial return, the conditions of this
 * Licence are imposed upon the receiver of the product, and all original and amended source code is included in any
 * transmitted product. You may be held legally responsible for any copyright infringement that is caused or encouraged by
 * your failure to abide by these terms and conditions.
 *
 * You are not permitted under this Licence to use this Software commercially. Use for which any financial return is received
 * shall be defined as commercial use, and includes:
 *
 * (1) integration of all or part of the source code or the Software into a product for sale or license by or on behalf of Licensee to third parties or
 *
 * (2) use of the Software or any derivative of it for research with the final aim of developing software products for sale or license to a third party or
 *
 * (3) use of the Software or any derivative of it for research with the final aim of developing non-software products for sale or license to a third party, or
 *
 * (4) use of the Software to provide any service to an external organisation for which payment is received.
 *
 * If you are interested in using the Software commercially, please contact Isis Innovation Limited ("Isis"), the technology transfer company of the University, to negotiate a licence.
 *
 * Contact details are: innovation@isis.ox.ac.uk quoting reference MG/11431.
 * */

//NikitaGourianov: Need to include a few C++ standard libraries Dr. Lubasch's PEPS library makes use of.
#include <cstdlib>
#include <typeinfo>
#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <string>
#include <sstream>

//NikitaGourianov: Need to connect with the TNT Library to make use of its tier 0 routines.
#ifdef __cplusplus
extern "C" {
#endif
    
#include "tnt.h"
#include "../lib/headers/tnt_int.h"
#include "../lib/headers/dec_la.h"
    
#ifdef __cplusplus
}
#endif

//NikitaGourianov: Now include Dr. Lubasch's (and modified by Nikita Gourianov) PEPS library
#include "../PEPSlib/MathAuxiliary.h"
#include "../PEPSlib/Vector.h"
#include "../PEPSlib/Tensor.h"
#include "../PEPSlib/Matrix.h"
#include "../PEPSlib/MPS.h"
#include "../PEPSlib/MPO.h"
#include "../PEPSlib/PEPS.h"
#include "../PEPSlib/PEPO.h"
#include "../PEPSlib/Operator2D.h"
#include "../PEPSlib/Hamiltonian2D.h"
#include "../PEPSlib/IsingModel2D.h"
#include "../PEPSlib/TimeEvolution2D.h"
