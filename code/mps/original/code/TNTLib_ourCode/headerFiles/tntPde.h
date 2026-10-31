/*
Authors: Nikita Gourianov
$LastChangedDate $
(c) University of Oxford 2020
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

/* include header for core library */
#include "tntMps.h"

/* QR decomposition implemented through SVD for testing purposes */
tntNode tntNodeQR_fake(tntNode A, tntLegLabel Qlegs, tntLegLabel legQR, tntLegLabel legRQ);

/* Supporting functions for solving PDEs */
void tntMpsCalcEntanSpec(tntNetwork psi, int chi, int orth_centre, tntDoubleArray *entanSpec, tntDoubleArray *entanEntr);
void tntMpsPdeProlongate(tntNetwork* mps, int L_new, int chi, char BC[], double* errp);
double tntMpsProlongContract(tntNetwork mpsp, int chi);

tntNetwork tntMpsPde1DCreateOrder1Diff(char BC[], unsigned L);
tntNetwork tntMpsPde1DCreateOrder2Diff(char BC[], unsigned L);
tntNetwork tntMpsPde1DCreateOrder4Diff(char BC[], unsigned L);
tntNetwork tntMpsPde1DCreateKsOrder2PlusOrder4Operator(char BC[], unsigned L, double mu, double nu, double lengthScale);
tntNetwork tntMpsPde1DCreateProlongator(unsigned L, char BC[]);

tntNetwork tntMpsPde2DCreateProlongator(unsigned L, char BC[]);
tntNetwork tntMpsPde2DCreateOrder1Diff1(char BC[],unsigned L);
tntNetwork tntMpsPde2DCreateOrder1Diff2(char BC[],unsigned L);
tntNetwork tntMpsPde2DCreateOrder2Diff1(char BC[], unsigned L);
tntNetwork tntMpsPde2DCreateOrder2Diff2(char BC[], unsigned L);
tntNetwork tntMpsPde2DCreateCrossDiffs_d1xd1(char BC[], unsigned L);
tntNetwork tntMpsPde2DCreateCrossDiffs_d2xd2(char BC[], unsigned L);

tntNetwork tntMpsPde3DCreateProlongator(unsigned L, char BC[]); //UNTESTED.
tntNetwork tntMpsPde3DCreateOrder1Diff1(char BC[], unsigned L);
tntNetwork tntMpsPde3DCreateOrder1Diff2(char BC[], unsigned L);
tntNetwork tntMpsPde3DCreateOrder1Diff3(char BC[], unsigned L);
tntNetwork tntMpsPde3DCreateOrder2Diff1(char BC[], unsigned L);
tntNetwork tntMpsPde3DCreateOrder2Diff2(char BC[], unsigned L);
tntNetwork tntMpsPde3DCreateOrder2Diff3(char BC[], unsigned L);
tntNetwork tntMpsPde3DCreateCrossDiffs_d1xd1(char BC[], unsigned L);
tntNetwork tntMpsPde3DCreateCrossDiffs_d2xd2(char BC[], unsigned L);
tntNetwork tntMpsPde3DCreateCrossDiffs_d3xd3(char BC[], unsigned L);

void tntMpsPdeMpsInit(tntNetwork u, tntNetwork v, tntNodeArray *betas, tntNodeArray *gammas);
void tntMpsPdeMpoInit(tntNetwork u, tntNetwork o, tntNetwork v, tntNodeArray *betas, tntNodeArray *gammas);
void tntMpsPdeMpsUpdate(tntNode U, tntNode V, tntLegLabel dir, unsigned j, tntNodeArray *Oeff);
void tntMpsPdeMpoUpdate(tntNode U, tntNode O, tntNode V, tntLegLabel dir, unsigned j, tntNodeArray *Oeff);

double tntMpoTruncate(tntNetwork mpo,int chi);
tntNetwork tntMps1MpoMps2Connect(tntNetwork mps1, tntNetwork mpo, tntNetwork mps2);
tntComplex tntMps1MpoMps2Product(tntNetwork mps1, tntNetwork mpo, tntNetwork mps2);
tntNetwork tntMpsHadamard(tntNetwork mps1, tntNetwork mps2, double epsilon);
tntNetwork tntMpsCreateHadamardMpo(tntNetwork mps);
void tntMpsHadamard1sFullSweep(tntNetwork hadProd, tntNetwork mps1, tntNetwork mps2I, tntNodeArray* Rxmps2I_left, tntNodeArray* Rxmps2I_right);
void tntMpsHadamard1sSweep(tntNetwork hadProd, tntNetwork mps1,tntNetwork mps2I,tntNodeArray *Rxmps2I_left, tntNodeArray *Rxmps2I_right,
                           tntLegLabel dir, int start_site, int end_site);
void tntHadamard1sLocalUpdate(tntNode *A, tntNode B, tntNode O, tntNode Rxmps2I_left, tntNode Rxmps2I_righ);
tntNode tntNodePseudoInverse(tntNode A, tntLegLabel inLegs, tntLegLabel outLegs, tntLegLabel legMap, double tol);
void tntNetworkArrayFree(tntNetwork *networkArray, unsigned arrayLength);
tntNode tntNodeAntiDirectSum(tntNode A, tntNode B, tntLegLabel expandedLegs);
tntNode tntNetworkPropagateSite_extLU(tntNetwork* HeffsumA, unsigned numHeffTerms, tntComplex dt, int dirTE);
tntNode tntNetworkPropagateSite_intLU(tntNode A, tntNode b);
