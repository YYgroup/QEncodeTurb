/*
Authors: Sarah Al-Assam, Stephen Clark, Dieter Jaksch (Oxford) and Chris Goodyer (NAG)
$LastChangedDate $
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

/* Function prototypes for all library functions available to users */

/* Core network level functions, listed in alphabetical order */

tntNetwork tntNetworkConjCopy(tntNetwork tnw);

tntNetwork tntNetworkCopy(tntNetwork tnw);

tntNetwork tntNetworkCreate(void);

void tntNetworkFree(tntNetwork *tnwp);

tntNode tntNetworkMinSite(tntNode tnA, tntNetwork *nwH, unsigned NH, tntNode (*eig_contract)(tntNode, tntNetwork), void (*eig_prep)(tntNode, tntNetwork), tntComplex *eigval);

tntNode tntNetworkPropagateSite(tntNetwork* HeffsumA, unsigned numHeffTerms, tntNode(*evolve_contract)(tntNode,tntNetwork), void (*evolve_prep)(tntNode* , tntNode*, tntNetwork), tntComplex dt, int dirTE, char method[]);

//void tntNetworkSchmidtVecInit(tntNetwork nwp, unsigned N, tntNode A);

tntNetwork tntNetworkSplit(tntNetwork tnw, 
                           tntNode tnFirst_1, tntLegLabel leg_start_1, tntNode tnLast_1, tntLegLabel leg_end_1,
                           tntNode tnFirst_2, tntLegLabel leg_start_2, tntNode tnLast_2, tntLegLabel leg_end_2,
                           unsigned numSplit_1, tntNode *tnSplit_1, unsigned numSplit_2, tntNode *tnSplit_2);

void tntNetworkToNodeGroup(tntNetwork *tnwp, int strip_connections);

/* Node functions, listed in alphabetical order */

void tntNodeAdd(tntNode tnA, tntNode tnB);

void tntNodeAddLeg(tntNode A, tntLegLabel legA);

void tntNodeClearQN(tntNode A);

/* Stringify for variable names */
#define tntNodeContract(tnA, ...) tntNodeContract_IMP(tnA, __VA_ARGS__, NULL, NULL, NULL)
#define tntNodeContract_IMP(tnA, tnB, legMapAC, legMapBC, ...) tntNodeNamedContract(#tnA, (tnA), #tnB, (tnB), (legMapAC), (legMapBC))
tntNode tntNodeNamedContract(const char *varnameA, tntNode tnA, const char *varnameB, tntNode tnB, tntLegLabel legMapAC, tntLegLabel legMapBC);

/* Fill up variable length list with TNT_MAXNUMNODES-2 NULL values to allow possibility of 2 to TNT_MAXNUMNODES nodes in contract list */
#define tntNodeListContract(remleglabels, ...) tntNodeListContract_IMP(remleglabels, __VA_ARGS__, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL) 
#define tntNodeContractList(remleglabels, ...) tntNodeListContract_IMP(remleglabels, __VA_ARGS__, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL) 
#define tntNodeListContract_IMP(remleglabels, tnA, tnB, tnC, tnD, tnE, tnF, tnG, tnH, tnI, tnJ, tnK, tnL, ...) tntNodeStaticListContract(remleglabels, #tnA, (tnA), #tnB, (tnB), #tnC, (tnC),  #tnD, (tnD), #tnE, (tnE), #tnF, (tnF),  #tnG, (tnG), #tnH, (tnH), #tnI, (tnI), #tnJ, (tnJ), #tnK, (tnK), #tnL, (tnL))
tntNode tntNodeStaticListContract(tntLegLabel remleglabels, 
                                  const char *varnameA, tntNode tnA, const char *varnameB, tntNode tnB, const char *varnameC, tntNode tnC, const char *varnameD, tntNode tnD, 
                                  const char *varnameE, tntNode tnE, const char *varnameF, tntNode tnF, const char *varnameG, tntNode tnG, const char *varnameH, tntNode tnH, 
                                  const char *varnameI, tntNode tnI, const char *varnameJ, tntNode tnJ, const char *varnameK, tntNode tnK, const char *varnameL, tntNode tnL);

#define tntNodeCopy(...) tntNodeCopy_IMP(__VA_ARGS__, 0, 0)
#define tntNodeCopy_IMP(tn,conj,...) tntNodeCopyConj(tn,conj)
tntNode tntNodeCopyConj(tntNode tn, int conj);

#define tntNodeCreate(...) tntNodeCreate_IMP(__VA_ARGS__, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#define tntNodeCreate_IMP(nodevals,leglabels,dima, dimb, dimc, dimd, dime, dimf, dimg, dimh, dimi, dimj, dimk, diml, ...) tntNodeCreateStatic(nodevals,leglabels,dima, dimb, dimc, dimd, dime, dimf, dimg, dimh, dimi, dimj, dimk, diml)
tntNode tntNodeCreateStatic(tntComplexArray *nodeVals, tntLegLabel leglabels, unsigned dima, unsigned dimb, unsigned dimc, unsigned dimd, unsigned dime, unsigned dimf, unsigned dimg, unsigned dimh, unsigned dimi, unsigned dimj, unsigned dimk, unsigned diml);

#define tntNodeFuncCreate(...) tntNodeFuncCreate_IMP(__VA_ARGS__, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#define tntNodeFuncCreate_IMP(nodevals,numops,opstr,leglabels,dima, dimb, dimc, dimd, dime, dimf, dimg, dimh, dimi, dimj, dimk, diml, ...) tntNodeFuncCreateStatic(nodevals,numops,opstr,leglabels,dima, dimb, dimc, dimd, dime, dimf, dimg, dimh, dimi, dimj, dimk, diml)
tntNode tntNodeFuncCreateStatic(tntComplexArray *nodeVals, unsigned numops, char *opstr, tntLegLabel leglabels, unsigned dima, unsigned dimb, unsigned dimc, unsigned dimd, unsigned dime, unsigned dimf, unsigned dimg, unsigned dimh, unsigned dimi, unsigned dimj, unsigned dimk, unsigned diml);

tntNode tntNodeCreateEyeLeg(tntNode tnA, tntLegLabel legA, tntNode tnB, tntLegLabel legB);

tntNode tntNodeCreateEyeOp(tntNode tn);

tntNode tntMpsCreateEye2(tntNode A, tntLegLabel outLegs, tntLegLabel inLegs);

tntNode tntNodeDirectSum(tntNode A, tntNode B, tntLegLabel leglabels);

void tntNodeExp(tntNode tn, tntComplex sc_val, tntLegLabel rowlegs, tntLegLabel collegs);

tntNode tntNodeFindConn(tntNode tn, tntLegLabel legA);

tntNode tntNodeFindFirst(tntNetwork tnw);

tntNode tntNodeFindLast(tntNetwork tnw);

void tntNodeFree(tntNode *tn);

void tntNodeReshape(tntNode *A, tntNode B);

void tntNodeReshapeAsMatrix(tntNode A, tntLegLabel rowlegs, tntLegLabel collegs);

unsigned tntNodeGetNumLegs(tntNode A);

void tntNodeFuse(tntNode A, tntLegLabel fuseLegs, tntLegLabel newLabelA, tntLegLabel newLabelB);

tntComplexArray tntNodeGetDiag(tntNode tn, tntLegLabel rowlegs, tntLegLabel collegs);

tntComplex tntNodeGetFirstEl(tntNode tn);

unsigned tntNodeGetLegDim(tntNode tn, tntLegLabel legid);

int tntNodeGetLegDir(tntNode tn, tntLegLabel legid);

tntComplexArray tntNodeGetMatrix(tntNode A, tntLegLabel rowlegs, tntLegLabel collegs);

double tntNodeGetNorm(tntNode tn);

tntIntArray tntNodeGetQN(tntNode tn, tntLegLabel legid);

tntComplex tntNodeGetTrace(tntNode tn, tntLegLabel rowlegs, tntLegLabel collegs);

void tntNodeGroupFree(tntNode *tn);

void tntNodeInsert(tntNode tnI, tntLegLabel legIA, tntLegLabel legIB, tntNode tnA, tntLegLabel legA, tntNode tnB, tntLegLabel legB);

void tntNodeInsertAtStart(tntNode tnI, tntLegLabel legInwstart, tntLegLabel legIfirst, tntNetwork tnw);

void tntNodeInsertAtEnd(tntNode tnI, tntLegLabel legIlast, tntLegLabel legInwend, tntNetwork tnw);

void tntNodeInverseFact(tntNode tn, tntLegLabel legid);

unsigned tntNodeIsCovariant(tntNode tn);

unsigned tntNodeIsFunctional(tntNode tn);

void tntNodeJoin(struct tnode *tnA, tntLegLabel legidA, struct tnode *tnB, tntLegLabel legidB);

void tntNodeMapLegs(tntNode tn, tntLegLabel legmap);

void tntNodeMakeCovariantQN(tntNode tn);

void tntNodePrepContract(tntNode tnA, tntNode tnB);

/* Stringify for variable names */
#define tntNodePrintInfo(name) tntNodePrintNamedInfo(#name, (name), 0)
#define tntNodePrintAll(name) tntNodePrintNamedInfo(#name, (name), 1)
void tntNodePrintNamedInfo(const char *varname, tntNode tn, int printvalues);

#define tntNodePrintAsMatrix(name, legsA,  legsB) tntNodePrintNamedAsMatrix(#name, (name), legsA, legsB)
void tntNodePrintNamedAsMatrix(const char *varname, tntNode tn, tntLegLabel rowlegs, tntLegLabel collegs);

#define tntNodeRemove(...) tntNodeRemove_IMP(__VA_ARGS__,NULL)
#define tntNodeRemove_IMP(A, connLegs, ...) tntNodeRemove_(A, connLegs)
void tntNodeRemove_(tntNode A, tntLegLabel connLegs);

void tntNodeReplace(tntNode tnA, tntNode tnB);

void tntNodeScaleReal(tntNode tn, double val);

void tntNodeScaleComplex(tntNode tn, tntComplex val);

void tntNodeSetCompParam(tntNode tn, tntComplex pval, unsigned pnum);

/* Stringify for variable names */
#define tntNodeContractSelf(tn) tntNodeNamedSelfContract(#tn, (tn))
#define tntNodeSelfContract(tn) tntNodeNamedSelfContract(#tn, (tn))
void tntNodeNamedSelfContract(const char *varname, tntNode *tn);

void tntNodeSetRealParam(tntNode tn, double pval, unsigned pnum);

void tntNodeSetQN(tntNode tn, tntLegLabel legid, tntIntArray *qvals, int legdir);

void tntNodeSetQN_nowarn(tntNode tn, tntLegLabel legid, tntIntArray *qvals, int legdir);

void tntNodeSplit(tntNode tnA, tntNode tnB);

void tntNodeSqueeze(tntNode tn, tntLegLabel squeezeLegs);

#define tntNodeSVD(...) tntNodeSVD_IMP(__VA_ARGS__,NULL,NULL)
#define tntNodeSVD_IMP(A, Ulegs, legUS, legSU, legSV, legVS, chi, err, newlegs, ...) tntNodeSVD_((A), #A, Ulegs, legUS, legSU, legSV, legVS, chi, err, newlegs)
tntNode tntNodeSVD_(tntNode A, const char *varnameA, tntLegLabel Ulegs, tntLegLabel legUS, tntLegLabel legSU, tntLegLabel legSV, tntLegLabel legVS, int chi, double * err, tntLegLabel legmap);

#define tntNodeQR(...) tntNodeQR_IMP(__VA_ARGS__,NULL)
#define tntNodeQR_IMP(A, Qlegs, legQR, legRQ, newlegs, ...) tntNodeQR_((A), #A, Qlegs, legQR, legRQ, newlegs)
tntNode tntNodeQR_(tntNode A, const char *varnameA, tntLegLabel Qlegs, tntLegLabel legQR, tntLegLabel legRQ, tntLegLabel legmap);

/* New load functions */
/* Expand out to the full number of arguments */
#define tntIntParamsLoad(fname, ...) tntIntParamsLoad_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntIntParamsNamedLoad(fname, ...) tntIntParamsNamedLoad_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntIntParamsLoad_IMP(fname, A, B, C, D, E, ...) tntIntParamsLoadStatic(fname, &(A), #A, &(B), #B, &(C), #C, &(D), #D, &(E), #E)
#define tntIntParamsNamedLoad_IMP(fname, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntIntParamsLoadStatic(fname, &(A), nA, &(B), nB, &(C), nC, &(D), nD, &(E), nE)
void tntIntParamsLoadStatic(const char *fname, int *varA, const char *varnameA, int *varB, const char *varnameB, 
                            int *varC, const char *varnameC, int *varD, const char *varnameD, int *varE, const char *varnameE);

/* Expand out to the full number of arguments */
#define tntDoubleParamsLoad(fname, ...) tntDoubleParamsLoad_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntDoubleParamsNamedLoad(fname, ...) tntDoubleParamsNamedLoad_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntDoubleParamsLoad_IMP(fname, A, B, C, D, E, ...) tntDoubleParamsLoadStatic(fname, &(A), #A, &(B), #B, &(C), #C, &(D), #D, &(E), #E)
#define tntDoubleParamsNamedLoad_IMP(fname, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntDoubleParamsLoadStatic(fname, &(A), nA, &(B), nB, &(C), nC, &(D), nD, &(E), nE)
void tntDoubleParamsLoadStatic(const char *fname, double *varA, const char *varnameA, double *varB, const char *varnameB, 
                               double *varC, const char *varnameC, double *varD, const char *varnameD, double *varE, const char *varnameE);

/* Expand out to the full number of arguments */
#define tntComplexParamsLoad(fname, ...) tntComplexParamsLoad_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntComplexParamsNamedLoad(fname, ...) tntComplexParamsNamedLoad_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */ 
#define tntComplexParamsLoad_IMP(fname, A, B, C, D, E, ...) tntComplexParamsLoadStatic(fname, &(A), #A, &(B), #B, &(C), #C, &(D), #D, &(E), #E)
#define tntComplexParamsNamedLoad_IMP(fname, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntComplexParamsLoadStatic(fname, &(A), nA, &(B), nB, &(C), nC, &(D), nD, &(E), nE)
void tntComplexParamsLoadStatic(const char *fname, tntComplex *varA, const char *varnameA, tntComplex *varB, const char *varnameB, 
                                tntComplex *varC, const char *varnameC, tntComplex *varD, const char *varnameD, tntComplex *varE, const char *varnameE);

/* Expand out to the full number of arguments */
#define tntStringsLoad(fname, ...) tntStringsLoad_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntStringsNamedLoad(fname, ...) tntStringsNamedLoad_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntStringsLoad_IMP(fname, A, B, C, D, E, ...) tntStringsLoadStatic(fname, (A), #A, (B), #B, (C), #C, (D), #D, (E), #E)
#define tntStringsNamedLoad_IMP(fname, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntStringsLoadStatic(fname, (A), nA, (B), nB, (C), nC, (D), nD, (E), nE)
void tntStringsLoadStatic(const char *fname, char *varA, const char *varnameA, char *varB, const char *varnameB, 
                          char *varC, const char *varnameC, char *varD, const char *varnameD, char *varE, const char *varnameE);

/* Expand out to the full number of arguments */
#define tntIntArraysLoad(fname, ...) tntIntArraysLoad_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntIntArraysNamedLoad(fname, ...) tntIntArraysNamedLoad_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntIntArraysLoad_IMP(fname, A, B, C, D, E, ...) tntIntArraysLoadStatic(fname, &(A), #A, &(B), #B, &(C), #C, &(D), #D, &(E), #E)
#define tntIntArraysNamedLoad_IMP(fname, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntIntArraysLoadStatic(fname, &(A), nA, &(B), nB, &(C), nC, &(D), nD, &(E), nE)
void tntIntArraysLoadStatic(const char *fname, tntIntArray *arrA, const char *varnameA, tntIntArray *arrB, const char *varnameB, 
                            tntIntArray *arrC, const char *varnameC, tntIntArray *arrD, const char *varnameD, tntIntArray *arrE, const char *varnameE);
/* Expand out to the full number of arguments */
#define tntDoubleArraysLoad(fname, ...) tntDoubleArraysLoad_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntDoubleArraysNamedLoad(fname, ...) tntDoubleArraysNamedLoad_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntDoubleArraysLoad_IMP(fname, A, B, C, D, E, ...) tntDoubleArraysLoadStatic(fname, &(A), #A, &(B), #B, &(C), #C, &(D), #D, &(E), #E)
#define tntDoubleArraysNamedLoad_IMP(fname, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntDoubleArraysLoadStatic(fname, &(A), nA, &(B), nB, &(C), nC, &(D), nD, &(E), nE)
void tntDoubleArraysLoadStatic(const char *fname, tntDoubleArray *arrA, const char *varnameA, tntDoubleArray *arrB, const char *varnameB, 
                               tntDoubleArray *arrC, const char *varnameC, tntDoubleArray *arrD, const char *varnameD, tntDoubleArray *arrE, const char *varnameE);

/* Expand out to the full number of arguments */
#define tntComplexArraysLoad(fname, ...) tntComplexArraysLoad_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntComplexArraysNamedLoad(fname, ...) tntComplexArraysNamedLoad_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntComplexArraysLoad_IMP(fname, A, B, C, D, E, ...) tntComplexArraysLoadStatic(fname, &(A), #A, &(B), #B, &(C), #C, &(D), #D, &(E), #E)
#define tntComplexArraysNamedLoad_IMP(fname, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntComplexArraysLoadStatic(fname, &(A), nA, &(B), nB, &(C), nC, &(D), nD, &(E), nE)
void tntComplexArraysLoadStatic(const char *fname, tntComplexArray *arrA, const char *varnameA, tntComplexArray *arrB, const char *varnameB, 
                                tntComplexArray *arrC, const char *varnameC, tntComplexArray *arrD, const char *varnameD, tntComplexArray *arrE, const char *varnameE);

/* Expand out to the full number of arguments */
#define tntNodesLoad(fname, ...) tntNodesLoad_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntNodesNamedLoad(fname, ...) tntNodesNamedLoad_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntNodesLoad_IMP(fname, A, B, C, D, E, ...) tntNodesLoadStatic(fname, &(A), #A, &(B), #B, &(C), #C, &(D), #D, &(E), #E)
#define tntNodesNamedLoad_IMP(fname, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntNodesLoadStatic(fname, &(A), nA, &(B), nB, &(C), nC, &(D), nD, &(E), nE)
void tntNodesLoadStatic(const char *fname, tntNode *varA, const char *varnameA, tntNode *varB, const char *varnameB, 
                        tntNode *varC, const char *varnameC, tntNode *varD, const char *varnameD, tntNode *varE, const char *varnameE);

/* Expand out to the full number of arguments */
#define tntNodeArraysLoad(fname, ...) tntNodeArraysLoad_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntNodeArraysNamedLoad(fname, ...) tntNodeArraysNamedLoad_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntNodeArraysLoad_IMP(fname, A, B, C, D, E, ...) tntNodeArraysLoadStatic(fname, &(A), #A, &(B), #B, &(C), #C, &(D), #D, &(E), #E)
#define tntNodeArraysNamedLoad_IMP(fname, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntNodeArraysLoadStatic(fname, &(A), nA, &(B), nB, &(C), nC, &(D), nD, &(E), nE)
void tntNodeArraysLoadStatic(const char *fname, tntNodeArray *arrA, const char *varnameA, tntNodeArray *arrB, const char *varnameB, 
                             tntNodeArray *arrC, const char *varnameC, tntNodeArray *arrD, const char *varnameD, tntNodeArray *arrE, const char *varnameE);

/* Expand out to the full number of arguments */
#define tntExOpsLoad(fname, ...) tntExOpsLoad_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntExOpsNamedLoad(fname, ...) tntExOpsNamedLoad_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntExOpsLoad_IMP(fname, A, B, C, D, E, ...) tntExOpsLoadStatic(fname, &(A), #A, &(B), #B, &(C), #C, &(D), #D, &(E), #E)
#define tntExOpsNamedLoad_IMP(fname, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntExOpsLoadStatic(fname, &(A), nA, &(B), nB, &(C), nC, &(D), nD, &(E), nE)
void tntExOpsLoadStatic(const char *fname, tntExOp *varA, const char *varnameA, tntExOp *varB, const char *varnameB, 
                        tntExOp *varC, const char *varnameC, tntExOp *varD, const char *varnameD, tntExOp *varE, const char *varnameE);

/* Expand out to the full number of arguments */
#define tntNetworksLoad(fname, ...) tntNetworksLoad_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntNetworksNamedLoad(fname, ...) tntNetworksNamedLoad_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntNetworksLoad_IMP(fname, A, B, C, D, E, ...) tntNetworksLoadStatic(fname, &(A), #A, &(B), #B, &(C), #C, &(D), #D, &(E), #E)
#define tntNetworksNamedLoad_IMP(fname, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntNetworksLoadStatic(fname, &(A), nA, &(B), nB, &(C), nC, &(D), nD, &(E), nE)
void tntNetworksLoadStatic(const char *fname, tntNetwork *varA, const char *varnameA, tntNetwork *varB, const char *varnameB, 
                           tntNetwork *varC, const char *varnameC, tntNetwork *varD, const char *varnameD, tntNetwork *varE, const char *varnameE);

/* New save functions */
/* Expand out to the full number of arguments */
#define tntIntParamsSave(fname, ...) tntIntParamsSave_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntIntParamsNamedSave(fname, ...) tntIntParamsNamedSave_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntIntParamsSave_IMP(fname, A, B, C, D, E, ...) tntIntParamsSaveStatic(fname, A, #A, B, #B, C, #C, D, #D, E, #E)
#define tntIntParamsNamedSave_IMP(fname, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntIntParamsSaveStatic(fname, A, nA, B, nB, C, nC, D, nD, E, nE)
void tntIntParamsSaveStatic(const char *fname, int varA, const char *varnameA, int varB, const char *varnameB, 
                            int varC, const char *varnameC, int varD, const char *varnameD, int varE, const char *varnameE);

/* Expand out to the full number of arguments */
#define tntDoubleParamsSave(fname, ...) tntDoubleParamsSave_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntDoubleParamsNamedSave(fname, ...) tntDoubleParamsNamedSave_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntDoubleParamsSave_IMP(fname, A, B, C, D, E, ...) tntDoubleParamsSaveStatic(fname, A, #A, B, #B, C, #C, D, #D, E, #E)
#define tntDoubleParamsNamedSave_IMP(fname, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntDoubleParamsSaveStatic(fname, A, nA, B, nB, C, nC, D, nD, E, nE)
void tntDoubleParamsSaveStatic(const char *fname, double varA, const char *varnameA, double varB, const char *varnameB, 
                               double varC, const char *varnameC, double varD, const char *varnameD, double varE, const char *varnameE);

/* Expand out to the full number of arguments */
#define tntComplexParamsSave(fname, ...) tntComplexParamsSave_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntComplexParamsNamedSave(fname, ...) tntComplexParamsNamedSave_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntComplexParamsSave_IMP(fname, A, B, C, D, E, ...) tntComplexParamsSaveStatic(fname, A, #A, B, #B, C, #C, D, #D, E, #E)
#define tntComplexParamsNamedSave_IMP(fname, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntComplexParamsSaveStatic(fname, A, nA, B, nB, C, nC, D, nD, E, nE)
void tntComplexParamsSaveStatic(const char *fname, tntComplex varA, const char *varnameA, tntComplex varB, const char *varnameB, 
                                tntComplex varC, const char *varnameC, tntComplex varD, const char *varnameD, tntComplex varE, const char *varnameE);

/* Expand out to the full number of arguments */
#define tntStringsSave(fname, ...) tntStringsSave_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntStringsNamedSave(fname, ...) tntStringsNamedSave_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntStringsSave_IMP(fname, A, B, C, D, E, ...) tntStringsSaveStatic(fname, A, #A, B, #B, C, #C, D, #D, E, #E)
#define tntStringsNamedSave_IMP(fname, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntStringsSaveStatic(fname, A, nA, B, nB, C, nC, D, nD, E, nE)
void tntStringsSaveStatic(const char *fname, const char *varA, const char *varnameA, const char *varB, const char *varnameB, 
                          const char *varC, const char *varnameC, const char *varD, const char *varnameD, const char *varE, const char *varnameE);

/* Expand out to the full number of arguments */
#define tntIntArraysSave(fname, ...) tntIntArraysSave_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntIntArraysNamedSave(fname, ...) tntIntArraysNamedSave_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntIntArraysSave_IMP(fname, A, B, C, D, E, ...) tntIntArraysSaveStatic(fname, A, #A, B, #B, C, #C, D, #D, E, #E)
#define tntIntArraysNamedSave_IMP(fname, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntIntArraysSaveStatic(fname, A, nA, B, nB, C, nC, D, nD, E, nE)
void tntIntArraysSaveStatic(const char *fname, tntIntArray arrA, const char *varnameA, tntIntArray arrB, const char *varnameB, 
                            tntIntArray arrC, const char *varnameC, tntIntArray arrD, const char *varnameD, tntIntArray arrE, const char *varnameE);

/* Expand out to the full number of arguments */
#define tntDoubleArraysSave(fname, ...) tntDoubleArraysSave_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntDoubleArraysNamedSave(fname, ...) tntDoubleArraysNamedSave_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntDoubleArraysSave_IMP(fname, A, B, C, D, E, ...) tntDoubleArraysSaveStatic(fname, A, #A, B, #B, C, #C, D, #D, E, #E)
#define tntDoubleArraysNamedSave_IMP(fname, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntDoubleArraysSaveStatic(fname, A, nA, B, nB, C, nC, D, nD, E, nE)
void tntDoubleArraysSaveStatic(const char *fname, tntDoubleArray arrA, const char *varnameA, tntDoubleArray arrB, const char *varnameB, 
                               tntDoubleArray arrC, const char *varnameC, tntDoubleArray arrD, const char *varnameD, tntDoubleArray arrE, const char *varnameE);

/* Expand out to the full number of arguments */
#define tntComplexArraysSave(fname, ...) tntComplexArraysSave_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntComplexArraysNamedSave(fname, ...) tntComplexArraysNamedSave_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntComplexArraysSave_IMP(fname, A, B, C, D, E, ...) tntComplexArraysSaveStatic(fname, A, #A, B, #B, C, #C, D, #D, E, #E)
#define tntComplexArraysNamedSave_IMP(fname, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntComplexArraysSaveStatic(fname, A, nA, B, nB, C, nC, D, nD, E, nE)
void tntComplexArraysSaveStatic(const char *fname, tntComplexArray arrA, const char *varnameA, tntComplexArray arrB, const char *varnameB, 
                                tntComplexArray arrC, const char *varnameC, tntComplexArray arrD, const char *varnameD, tntComplexArray arrE, const char *varnameE);

/* Expand out to the full number of arguments */
#define tntNodesSave(fname, ...) tntNodesSave_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntNodesNamedSave(fname, ...) tntNodesNamedSave_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntNodesSave_IMP(fname, A, B, C, D, E, ...) tntNodesSaveStatic(fname, A, #A, B, #B, C, #C, D, #D, E, #E)
#define tntNodesNamedSave_IMP(fname, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntNodesSaveStatic(fname, A, nA, B, nB, C, nC, D, nD, E, nE)
void tntNodesSaveStatic(const char *fname, tntNode varA, const char *varnameA, tntNode varB, const char *varnameB, 
                        tntNode varC, const char *varnameC, tntNode varD, const char *varnameD, tntNode varE, const char *varnameE);

/* Expand out to the full number of arguments */
#define tntNodeArraysSave(fname, ...) tntNodeArraysSave_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntNodeArraysNamedSave(fname, ...) tntNodeArraysNamedSave_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntNodeArraysSave_IMP(fname, A, B, C, D, E, ...) tntNodeArraysSaveStatic(fname, A, #A, B, #B, C, #C, D, #D, E, #E)
#define tntNodeArraysNamedSave_IMP(fname, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntNodeArraysSaveStatic(fname, A, nA, B, nB, C, nC, D, nD, E, nE)
void tntNodeArraysSaveStatic(const char *fname, tntNodeArray arrA, const char *varnameA, tntNodeArray arrB, const char *varnameB, 
                             tntNodeArray arrC, const char *varnameC, tntNodeArray arrD, const char *varnameD, tntNodeArray arrE, const char *varnameE);

/* Expand out to the full number of arguments */
#define tntExOpsSave(fname, ...) tntExOpsSave_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntExOpsNamedSave(fname, ...) tntExOpsNamedSave_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntExOpsSave_IMP(fname, A, B, C, D, E, ...) tntExOpsSaveStatic(fname, A, #A, B, #B, C, #C, D, #D, E, #E)
#define tntExOpsNamedSave_IMP(fname, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntExOpsSaveStatic(fname, A, nA, B, nB, C, nC, D, nD, E, nE)
void tntExOpsSaveStatic(const char *fname, tntExOp varA, const char *varnameA, tntExOp varB, const char *varnameB, 
                        tntExOp varC, const char *varnameC, tntExOp varD, const char *varnameD, tntExOp varE, const char *varnameE);

/* Expand out to the full number of arguments */
#define tntNetworksSave(fname, ...) tntNetworksSave_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntNetworksNamedSave(fname, ...) tntNetworksNamedSave_IMP(fname, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntNetworksSave_IMP(fname, A, B, C, D, E, ...) tntNetworksSaveStatic(fname, A, #A, B, #B, C, #C, D, #D, E, #E)
#define tntNetworksNamedSave_IMP(fname, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntNetworksSaveStatic(fname, A, nA, B, nB, C, nC, D, nD, E, nE)
void tntNetworksSaveStatic(const char *fname, tntNetwork varA, const char *varnameA, tntNetwork varB, const char *varnameB, 
                           tntNetwork varC, const char *varnameC, tntNetwork varD, const char *varnameD, tntNetwork varE, const char *varnameE);

/* New update functions */
/* Expand out to the full number of arguments */
#define tntIntParamsUpdate(fname, posn, ...) tntIntParamsUpdate_IMP(fname, posn, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntIntParamsNamedUpdate(fname, posn, ...) tntIntParamsNamedUpdate_IMP(fname, posn, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntIntParamsUpdate_IMP(fname, posn, A, B, C, D, E, ...) tntIntParamsUpdateStatic(fname, posn, A, #A, B, #B, C, #C, D, #D, E, #E)
#define tntIntParamsNamedUpdate_IMP(fname, posn, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntIntParamsUpdateStatic(fname, posn, A, nA, B, nB, C, nC, D, nD, E, nE)
void tntIntParamsUpdateStatic(const char *fname, unsigned posn, int varA, const char *varnameA, int varB, const char *varnameB, 
                              int varC, const char *varnameC, int varD, const char *varnameD, int varE, const char *varnameE);

/* Expand out to the full number of arguments */
#define tntDoubleParamsUpdate(fname, posn, ...) tntDoubleParamsUpdate_IMP(fname, posn, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntDoubleParamsNamedUpdate(fname, posn, ...) tntDoubleParamsNamedUpdate_IMP(fname, posn, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntDoubleParamsUpdate_IMP(fname, posn, A, B, C, D, E, ...) tntDoubleParamsUpdateStatic(fname, posn, A, #A, B, #B, C, #C, D, #D, E, #E)
#define tntDoubleParamsNamedUpdate_IMP(fname, posn, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntDoubleParamsUpdateStatic(fname, posn, A, nA, B, nB, C, nC, D, nD, E, nE)
void tntDoubleParamsUpdateStatic(const char *fname, unsigned posn, double varA, const char *varnameA, double varB, const char *varnameB, 
                                 double varC, const char *varnameC, double varD, const char *varnameD, double varE, const char *varnameE);

/* Expand out to the full number of arguments */
#define tntComplexParamsUpdate(fname, posn, ...) tntComplexParamsUpdate_IMP(fname, posn, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntComplexParamsNamedUpdate(fname, posn, ...) tntComplexParamsNamedUpdate_IMP(fname, posn, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntComplexParamsUpdate_IMP(fname, posn, A, B, C, D, E, ...) tntComplexParamsUpdateStatic(fname, posn, A, #A, B, #B, C, #C, D, #D, E, #E)
#define tntComplexParamsNamedUpdate_IMP(fname, posn, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntComplexParamsUpdateStatic(fname, posn, A, nA, B, nB, C, nC, D, nD, E, nE)
void tntComplexParamsUpdateStatic(const char *fname, unsigned posn, tntComplex varA, const char *varnameA, tntComplex varB, const char *varnameB, 
                                  tntComplex varC, const char *varnameC, tntComplex varD, const char *varnameD, tntComplex varE, const char *varnameE);

/* Expand out to the full number of arguments */
#define tntIntArraysUpdate(fname, posn, ...) tntIntArraysUpdate_IMP(fname, posn, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntIntArraysNamedUpdate(fname, posn, ...) tntIntArraysNamedUpdate_IMP(fname, posn, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntIntArraysUpdate_IMP(fname, posn, A, B, C, D, E, ...) tntIntArraysUpdateStatic(fname, posn, A, #A, B, #B, C, #C, D, #D, E, #E)
#define tntIntArraysNamedUpdate_IMP(fname, posn, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntIntArraysUpdateStatic(fname, posn, A, nA, B, nB, C, nC, D, nD, E, nE)
void tntIntArraysUpdateStatic(const char *fname, unsigned posn, tntIntArray arrA, const char *varnameA, tntIntArray arrB, const char *varnameB, 
                              tntIntArray arrC, const char *varnameC, tntIntArray arrD, const char *varnameD, tntIntArray arrE, const char *varnameE);

/* Expand out to the full number of arguments */
#define tntDoubleArraysUpdate(fname, posn, ...) tntDoubleArraysUpdate_IMP(fname, posn, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntDoubleArraysNamedUpdate(fname, posn, ...) tntDoubleArraysNamedUpdate_IMP(fname, posn, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntDoubleArraysUpdate_IMP(fname, posn, A, B, C, D, E, ...) tntDoubleArraysUpdateStatic(fname, posn, A, #A, B, #B, C, #C, D, #D, E, #E)
#define tntDoubleArraysNamedUpdate_IMP(fname, posn, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntDoubleArraysUpdateStatic(fname, posn, A, nA, B, nB, C, nC, D, nD, E, nE)
void tntDoubleArraysUpdateStatic(const char *fname, unsigned posn, tntDoubleArray arrA, const char *varnameA, tntDoubleArray arrB, const char *varnameB, 
                                 tntDoubleArray arrC, const char *varnameC, tntDoubleArray arrD, const char *varnameD, tntDoubleArray arrE, const char *varnameE);

/* Expand out to the full number of arguments */
#define tntComplexArraysUpdate(fname, posn, ...) tntComplexArraysUpdate_IMP(fname, posn, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
#define tntComplexArraysNamedUpdate(fname, posn, ...) tntComplexArraysNamedUpdate_IMP(fname, posn, __VA_ARGS__, __VA_ARGS__ , __VA_ARGS__ , __VA_ARGS__, __VA_ARGS__) 
/* Put the arguments into the form required by the main function */
#define tntComplexArraysUpdate_IMP(fname, posn, A, B, C, D, E, ...) tntComplexArraysUpdateStatic(fname, posn, A, #A, B, #B, C, #C, D, #D, E, #E)
#define tntComplexArraysNamedUpdate_IMP(fname, posn, A, nA, B, nB, C, nC, D, nD, E, nE, ...) tntComplexArraysUpdateStatic(fname, posn, A, nA, B, nB, C, nC, D, nD, E, nE)
void tntComplexArraysUpdateStatic(const char *fname, unsigned posn, tntComplexArray arrA, const char *varnameA, tntComplexArray arrB, const char *varnameB, 
                                  tntComplexArray arrC, const char *varnameC, tntComplexArray arrD, const char *varnameD, tntComplexArray arrE, const char *varnameE);

/* Functions for altering or inspecting system environment */

void tntSysQNClearWarnOn(void);

void tntSysQNClearWarnOff(void);

void tntSysTypeSet(int sysnum);

int tntSysTypeGet(void);

void tntSysBasisOpSet(tntNode basisOp);

void tntSysBasisOpFree(void);

tntNode tntSysBasisOpGet(void);

void tntSymmTypeSet(const char *symmTypeStr, unsigned num_qn);

int tntSymmTypeGet(void);

void tntSymmTypeUnset(void);

double tntWalltimeGet(void);

void tntWalltimeSet(double t);

unsigned tntSymmNumGet(void);

/* Other functions */

void tntQuietModeOn(void);

void tntQuietModeOff(void);

#define tntTruncType tntSVDTruncType
void tntSVDTruncType(const char *funcname);

void tntMaxEigIterSet(int maxnum);

void tntPrintVersion(void);

void tntSVDRoutineSet(int svdtype);

void tntlibraryPreference_CPUvsGPUSet(int libPref);

double tntSVDTolGet(void);

void tntSVDTolSet(double inpTOL);

double tntSVDTruncTolGet(void);

void tntSVDTruncTolSet(double truncTol);

double tntSVDRelTruncTolGet(void);

void tntSVDRelTruncTolSet(double truncTol);

double tntSVDTruncErrTolGet(void);

void tntSVDTruncErrTolSet(double truncTol);

void tntReshapeReuseDisable(void);

void tntReshapeReuseEnable(void);

void tntClearReshapes(int chi);

void tntInitialize(void);

void tntFinalize(void);

void tntProcessCLOptions(int argc, char **argv, char *loadname, char *loadname_op, char *saveprefix, int *chi, int *outputState);

int tntPrintf(const char *format, ...);

void tntSysInfoPrint(void);

void tntWarningPrint(const char *msg, ...);

tntComplex tntAtoC(const char * str);

#define tntErrorPrint(...) tntNamedErrorPrint(__FUNCTION__, __VA_ARGS__)
void tntNamedErrorPrint(const char *funcname, const char *msg, ...);

double tntComplexToReal(tntComplex var);

/* Prototypes of array functions, listed in alphabetical order */

#define tntComplexArrayAlloc(...) tntComplexArrayAlloc_IMP(__VA_ARGS__, 1)
#define tntComplexArrayAlloc_IMP(numrows, numcols, ...) tntComplexMatrixAlloc(numrows, numcols)
tntComplexArray tntComplexMatrixAlloc(unsigned numrows, unsigned numcols);

void tntComplexArrayFree(tntComplexArray *arr);

/* Stringify for variable names */
#define tntComplexArrayPrint(name) tntNamedComplexArrayPrint(#name, &(name))
void tntNamedComplexArrayPrint(const char *varname, tntComplexArray *arr);

#define tntDoubleArrayAlloc(...) tntDoubleArrayAlloc_IMP(__VA_ARGS__, 1)
#define tntDoubleArrayAlloc_IMP(numrows, numcols, ...) tntDoubleMatrixAlloc(numrows, numcols)
tntDoubleArray tntDoubleMatrixAlloc(unsigned numrows, unsigned numcols);

void tntDoubleArrayFree(tntDoubleArray *arr);

/* Stringify for variable names */
#define tntDoubleArrayPrint(name) tntNamedDoubleArrayPrint(#name, &(name))
void tntNamedDoubleArrayPrint(const char *varname, tntDoubleArray *arr);

#define tntIntArrayAlloc(...) tntIntArrayAlloc_IMP(__VA_ARGS__, 1)
#define tntIntArrayAlloc_IMP(numrows, numcols, ...) tntIntMatrixAlloc(numrows, numcols)
tntIntArray tntIntMatrixAlloc(unsigned numrows, unsigned numcols);

void tntIntArrayFree(tntIntArray *arr);

/* Stringify for variable names */
#define tntIntArrayPrint(name) tntNamedIntArrayPrint(#name, &(name))
void tntNamedIntArrayPrint(const char *varname, tntIntArray *arr);

int tntComplexArrayIsReal(tntComplexArray *arr);

tntDoubleArray tntComplexArrayToReal(tntComplexArray *arr);

tntNodeArray tntNodeArrayAlloc(unsigned num_elems);

void tntNodeArrayFree(tntNodeArray *arr);

tntStringArray tntStringArrayAlloc(unsigned num_elems);

void tntStringArrayFree(tntStringArray *arr);

tntIntArray tntIntArrayCreate(const char *array_string);

tntDoubleArray tntDoubleArrayCreate(const char *array_string);

tntComplexArray tntComplexArrayCreate(const char *array_string);

tntIntArray tntIntArrayCreateDiag(const char *array_string, int offs);

tntDoubleArray tntDoubleArrayCreateDiag(const char *array_string, int offs);

tntComplexArray tntComplexArrayCreateDiag(const char *array_string, int offs);

#define tntNodeArrayCreate(...) tntNodeArrayCreate_IMP(__VA_ARGS__, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL) 
#define tntNodeArrayCreate_IMP(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,...) tntNodeArrayCreateStatic(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X)
tntNodeArray tntNodeArrayCreateStatic(tntNode A, tntNode B, tntNode C, tntNode D, tntNode E, tntNode F, tntNode G, tntNode H, tntNode I, tntNode J, tntNode K, tntNode L, tntNode M, tntNode N, tntNode O, tntNode P, tntNode Q, tntNode R, tntNode S, tntNode T, tntNode U, tntNode V, tntNode W, tntNode X);

#define tntStringArrayCreate(...) tntStringArrayCreate_IMP(__VA_ARGS__, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL) 
#define tntStringArrayCreate_IMP(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,...) tntStringArrayCreateStatic(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X)
tntStringArray tntStringArrayCreateStatic(const char *A, const char *B, const char *C, const char *D, const char *E, const char *F, const char *G, const char *H, const char *I, const char *J, const char *K, const char *L, const char *M, const char *N, const char *O, const char *P, const char *Q, const char *R, const char *S, const char *T, const char *U, const char *V, const char *W, const char *X);
