#ifndef TNT_MPS_DOT_H
#define TNT_MPS_DOT_H

/*
Authors: Sarah Al-Assam, Nikita Gourianov, Stephen Clark, Dieter Jaksch and Michael Lubasch (Oxford) and Chris Goodyer (NAG)
$LastChangedDate: 2018-09-25 11:52:49 +0100 (Tue, 25 Sep 2018) $
(c) University of Oxford 2014
*/

/* Software Licence 
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
#include "tnt.h"

/* Commonly used leg labels */
#define TNT_MPS_L "L"
#define TNT_MPS_R "R"
#define TNT_MPS_D "D"
#define TNT_MPS_U "U"
#define TNT_MPS_TUL "U"
#define TNT_MPS_TUR "V"
#define TNT_MPS_TDL "D"
#define TNT_MPS_TDR "E"

/* Code to denote the system is a spin system. */
#define TNT_MPS_SPIN 1

/* Code to denote the system is a boson system. */
#define TNT_MPS_BOSON 2

/* error code for any error generated in MPS public functions */
#define TNT_MPS_ERRCODE 2

/*! \ingroup mps 
 * A structure that includes information about all the expectation values that should be taken.
 */
#define tntMpsExOp tntExOp

void tntMpoExpecOutput(tntNetwork rho, tntMpsExOp *Op, unsigned printOutput, unsigned saveOutput, char *saveprefix, unsigned counter);  

tntNetwork tntMpoHC(tntNetwork O);

tntNetwork tntMpoMpoProduct(tntNetwork O, tntNetwork P);

tntComplex tntMpoMpoTrace(tntNetwork mpoO, tntNetwork mpoP);

void tntMpoPmpoProduct(tntNetwork mpo, tntNodeArray *op, tntIntArray *sitenum);

tntComplex tntMpoPmpoTrace(tntNetwork mpo, tntNodeArray *op, tntIntArray *sitenum);

void tntMpoPropST2scConnect(tntNetwork mpo, tntNetwork PropB, tntNetwork PropT);

double tntMpoPropSTscContract(tntNetwork mpoProp, int chi);

double tntMpoPropST2scProduct(tntNetwork mpo, tntNetwork PropB, tntNetwork PropT, int chi);

tntNode tntMpoPropTwoSite(tntNode *Mlp, tntNode *Mrp, int chi,  double *err);

tntComplex tntMpoTrace(tntNetwork mpo);

tntNode tntMpsVarMinMpo2sContract_alt(tntNode theta, tntNetwork nwMV);

tntNetwork tntMpsAdd(tntNetwork mpsA, tntNetwork mpsB, int orthresult);

void tntMpsCreateBosonOp(unsigned Nmax, tntNode *b, tntNode *bd, tntNode *n, tntNode *os_int, tntNode *eye);

tntNetwork tntMpsCreateConfig(unsigned L, const char *config);

#define tntMpsCreateMpo(...) tntMpsCreateMpo_IMP(__VA_ARGS__, 0)
#define tntMpsCreateMpo_IMP(L, nnL, nnR, nnparam, os, osparam, ignore_qn, ...) tntMpsCreateMpo_(L, nnL, nnR, nnparam, os, osparam, ignore_qn)
tntNetwork tntMpsCreateMpo_(unsigned L, tntNodeArray *nnL, tntNodeArray *nnR, tntComplexArray *nnparam, tntNodeArray *os, tntComplexArray *osparam, unsigned ignore_qn);

tntNetwork tntMpsCreateEyeMpo(unsigned L, tntNode basisOp);

tntNetwork tntMpsCreateProductMpo(unsigned L, tntNodeArray *op, tntIntArray *sitenum);

tntNodeArray tntMpsCreatePropArray(unsigned L, tntComplex dtc, tntNodeArray *nnL, tntNodeArray *nnR, tntComplexArray *nnparam, tntNodeArray *os, tntComplexArray *osparam);

tntNetwork tntMpsCreatePropST2sc(unsigned L, tntComplex dtc, tntNodeArray *nnL, tntNodeArray *nnR, tntComplexArray *nnparam, tntNodeArray *os, tntComplexArray *osparam);

void tntMpsCreateSpinOp(unsigned TwoS, tntNode *Sx, tntNode *Sy, tntNode *Sz, tntNode *Sp, tntNode *Sm, tntNode *eye);

tntNetwork tntMpsCreateRandom(unsigned L, unsigned D);

tntNetwork tntMpsCreateSymmRandom(unsigned L, int *qn);

tntNode tntMpsCreateTwoSiteOp(tntNodeArray *Lnodes, tntNodeArray *Rnodes, tntComplexArray *params);

void tntMpsExOpFree(tntMpsExOp *ExOp);

void tntMpsExpecOutput(tntNetwork wf, tntMpsExOp *Op, int orthCentre, unsigned printOutput, unsigned saveOutput, char *saveprefix, unsigned counter);

#define tntMpoLength(wf) tntMpsLength(wf)
unsigned tntMpsLength(tntNetwork wf);

tntIntArray tntMpsOpGetQN(tntNode op);

void tntMpsMpoConnect(tntNetwork mps, tntNetwork mpo);

double tntMpsMpoContract(tntNetwork mpsmpo, int chi);

double tntMpsMpoProduct(tntNetwork mps, tntNetwork mpo, int chi);

tntNetwork tntMpsMpoMpsConnect(tntNetwork mps, tntNetwork mpo);

tntComplex tntMpsMpoMpsContract(tntNetwork sandwich);

void tntMpsMpoMpsInitOrth(tntNetwork mps, tntNetwork mpo, tntNodeArray *beta, tntNodeArray *gammas);

void tntMpsMpoMpsInit(tntNetwork mps, tntNetwork mpo, tntNodeArray *beta, tntNodeArray *gammas);

tntComplex tntMpsMpoMpsPartProduct(tntNetwork mps, tntNetwork mpo, unsigned site_left, unsigned site_right, tntNode beta, tntNode gamma);

tntComplex tntMpsMpoMpsProduct(tntNetwork mps, tntNetwork mpo);

void tntMpsMpsConnect(tntNetwork mpsA, tntNetwork mpsB);


#define tntMpsMpsContract(...) tntMpsMpsContract_IMP(__VA_ARGS__, NULL, NULL)
#define tntMpsMpsContract_IMP(sandwich, site_left, site_right, betas, gammas, ...) tntMpsMpsContract_(sandwich, site_left, site_right, betas, gammas)
tntComplex tntMpsMpsContract_(tntNetwork sandwich, int site_left, int site_right, tntNode beta, tntNode gamma);

tntComplex tntMpsMpsProduct(tntNetwork mpsA, tntNetwork mpsB);

double tntMpsNodeNodeMpsProduct(tntNetwork wf, tntNode op, unsigned sitenum);

void tntMpsNodeProduct(tntNetwork wf, tntNode op, unsigned sitenum);

void tntMpsOrth(tntNetwork wf, unsigned orth_centre);

void tntMpsOrthNorm(tntNetwork mps,unsigned orth_centre); 

#define tntMpsPmpoMpsProduct(...) tntMpsPmpoMpsProduct_IMP(__VA_ARGS__, NULL, NULL)
#define tntMpsPmpoMpsProduct_IMP(mps, op, sitenum, orth_centre, betas, gammas, ...) tntMpsPmpoMpsProduct_(mps, op, sitenum, orth_centre, betas, gammas)
tntComplex tntMpsPmpoMpsProduct_(tntNetwork mps, tntNodeArray *op, tntIntArray *sitenum, int orth_centre, tntNodeArray *betas, tntNodeArray *gammas);

void tntMpsPmpoProduct(tntNetwork mps, tntNodeArray *op, tntIntArray *sitenum);

void tntMpsProcessDmrgCLOptions(int argc, char **argv, int *max_it, double *precision);

void tntMpsProcessExpecOptions(int argc, char **argv, tntMpsExOp *Op);

void tntMpsProcessSysOptions(int argc, char **argv, tntNetwork *mps, tntNetwork *mpo, tntNetwork *prop, tntMpsExOp *Op);

void tntMpsProcessTebdCLOptions(int argc, char **argv, int *numsteps, int *tbigstep, int *update_chi,  double *emax, int *delta_chi, int *max_chi);

void tntMpsProcessTdvpCLOptions(int argc, char **argv, int *numsteps, int *tbigstep, int *update_chi,  double *emax, int *delta_chi, int *max_chi, double *DtRe, double *DtIm);

tntNetwork tntMpsPropArrayToST2sc(tntNodeArray Proparr);

void tntMpsPropST2scConnect(tntNetwork mps, tntNetwork Prop);

double tntMpsPropSTscContract(tntNetwork mpsProp, int chi);

double tntMpsPropST2scProduct(tntNetwork mps, tntNetwork Prop, int chi);

tntNode tntMpsPropTwoSite(tntNode *Mlp, tntNode *Mrp, int chi,  double *err);

void tntMpsSelfInit(tntNetwork mps, tntNodeArray *betas, tntNodeArray *gammas);

double tntMpsSelfProduct(tntNetwork mps, int orth_centre);

tntNetwork tntMpsVarMinMpo2sBuild(tntNode beta, tntNode gamma, tntNode theta, tntNode tno);

tntNode tntMpsVarMinMpo2sContract(tntNode theta, tntNetwork nwMV);

tntNetwork tntMpsVarMinMpo1sBuild(tntNode beta, tntNode gamma, tntNode A, tntNode tno);

tntNetwork tntMpsVarMinMpo0sBuild(tntNode beta, tntNode gamma, tntNode C);

tntNode tntMpsVarMinMpo1sContract(tntNode A, tntNetwork nwMV);

tntNode tntMpsVarMinMpo0sContract(tntNode C, tntNetwork nwMV);

double tntMpsVarMinMpo1sFullContract(tntNode beta, tntNode gamma, tntNode A, tntNode O);

#define tntMpsVarMinMpoInit(mps, mpo, betas, gammas) tntMpsMpoMpsInit(mps, mpo, betas, gammas)

double tntMpsVarMinMpo2sSweep(tntNetwork mps, tntLegLabel direction, int chi, tntNetwork mpo, tntNodeArray *HeffL, tntNodeArray *HeffR, double *err, int start_site, int end_site);

double tntMpsVarMinMpo1sSweep(tntNetwork mps, tntLegLabel direction, int chi, tntNetwork mpo, tntNodeArray *HeffL, tntNodeArray *HeffR, double *err, int start_site, int end_site);

tntNode tntDmrg2sLocalOps(tntNode Al, tntNode Ar, tntNode Ol, tntNode Or, tntNode beta, tntNode gamma, int chi, double *errp, double *eigvp);

tntNode tntDmrg1sLocalOptimise(tntNode A, tntNode O, tntNode beta, tntNode gamma, double *eigvp);

tntNode tntDmrg1sSubspaceExpand(tntNode Al, tntNode Ar, tntNode Ol, tntNode Or, tntNode beta, tntNode gamma, tntLegLabel dir, int chi, double *errp);

#define tntMpsVarMinMpo2sStep(...) tntMpsVarMinMpo2sStep_IMP(__VA_ARGS__, 0)
#define tntMpsVarMinMpo2sStep_IMP(mps, chi, mpo, HeffL, HeffR, err, start_site, ...) tntMpsVarMinMpo2sStep_(mps, chi, mpo, HeffL, HeffR, err, start_site)
double tntMpsVarMinMpo2sStep_(tntNetwork mps, int chi, tntNetwork mpo, tntNodeArray *HeffL, tntNodeArray *HeffR, double *err, unsigned start_site);

#define tntMpsVarMinMpo1sStep(...) tntMpsVarMinMpo1sStep_IMP(__VA_ARGS__, 0)
#define tntMpsVarMinMpo1sStep_IMP(mps, chi, mpo, HeffL, HeffR, err, start_site, ...) tntMpsVarMinMpo1sStep_(mps, chi, mpo, HeffL, HeffR, err, start_site)
double tntMpsVarMinMpo1sStep_(tntNetwork mps, int chi, tntNetwork mpo, tntNodeArray *HeffL, tntNodeArray *HeffR, double *err, unsigned start_site);

void tntMpsVarMinMpoUpdate(tntNode tnA, tntNode tno, tntLegLabel dir, unsigned j, tntNodeArray *Heff);

/********************************* TDVP MPS FUNCTIONS ******************************************/
void tntMpsTdvp2sStep(tntNetwork mps, int chi, tntNetwork *mpo_sum, tntNodeArray *(HeffL_sum[]), tntNodeArray *(HeffR_sum[]), double *errp, tntComplex dt, int skipBackwardsTEStep, int numMpoTerms, char method[]);
void tntMpsTdvp1sStep(tntNetwork mps, tntNetwork *mpo_sum, tntNodeArray *(HeffL_sum[]), tntNodeArray *(HeffR_sum[]), tntComplex dt, int skipBackwardsTEStep, int numMpoTerms, char method[]);
void tntMpsTdvp2sSweep(tntNetwork mps, tntLegLabel dir, int chi, tntNetwork *mpo_sum, tntNodeArray *(HeffL_sum[]), tntNodeArray *(HeffR_sum[]), double *errp, int start_site, int end_site, tntComplex dt, int skipBackwardsTEStep, int numMpoTerms, char method[]);
void tntMpsTdvp1sSweep(tntNetwork mps, tntLegLabel dir, tntNetwork *mpo_sum, tntNodeArray *(HeffL_sum[]), tntNodeArray *(HeffR_sum[]), int start_site, int end_site, tntComplex dt, int skipBackwardsTEStep, int numMpoTerms, char method[]);
tntNode tntTdvp2sLocalTE(tntNode Al, tntNode Ar, tntNode *Ol_sum, tntNode *Or_sum, tntNode *beta_sum, tntNode *gamma_sum, int chi, double *errp, tntComplex dt, int dirTE, int numMpoTerms, char method[]);
tntNode tntTdvp1sLocalTE(tntNode A, tntNode *O_sum, tntNode *beta_sum, tntNode *gamma_sum, tntComplex dt, int dirTE, int numMpoTerms, char method[]);
tntNode tntTdvp0sLocalTE(tntNode C, tntNode *beta_sum, tntNode *gamma_sum, tntComplex dt, int dirTE, int numMpoTerms, char method[]);
tntNetwork tntMpsTdvp1sExpandInternalBasis(tntNetwork wfIn, int chi_max);
void tntMpsVarMinMpo2sPrepContract(tntNode *theta_contracted, tntNode* theta, tntNetwork Htheta);
void tntMpsVarMinMpo1sPrepContract(tntNode *A_contracted, tntNode *A, tntNetwork HA);
void tntMpsVarMinMpo0sPrepContract(tntNode *C_contracted, tntNode *C, tntNetwork KC);

/* Deprecated functions */
tntNodeArray tntMpsCreatePropagator(unsigned L, tntComplex dtc, tntNodeArray *nnL, tntNodeArray *nnR, tntComplexArray *nnparam, tntNodeArray *os, tntComplexArray *osparam);

double tntMpsExpecSing(tntNetwork wf, unsigned num_op, tntNode op[], unsigned sitenum[], int orth_centre, double norm);

void tntMpsFlip(tntNetwork mpsA, tntNetwork mpsB);

void tntMpsProcessInitOptions(int argc, char **argv, tntNetwork *wf, tntNetwork *H, tntNodeArray *Uarr, tntMpsExOp *Op);

double tntMpsOverlap(tntNetwork wfa, tntNetwork wfb);

double tntMpsMpoExpec(tntNetwork mps, tntNetwork mpo);

double tntMpsMpoProdExpec(tntNetwork wf, unsigned num_op, tntNode op[], unsigned sitenum[]);

double tntMpsNorm(tntNetwork wf, int orth_centre);

double tntMpsTruncate(tntNetwork mps, unsigned orth_centre, int chi);

double tntMpsTebdSweepInhom(tntNetwork wf, int direction, int chi, tntNodeArray U);

double tntMpsTebdSweepUni(tntNetwork wf, int direction, int chi, tntNode UL, tntNode UM, tntNode UR);

double tntMpsTebdSweepVar(tntNetwork wf, int direction, int chi, tntNode U, tntDoubleArray params);

#endif
