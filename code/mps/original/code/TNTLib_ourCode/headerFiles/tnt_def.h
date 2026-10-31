#ifndef TNT_DEFS_DOT_H
#define TNT_DEFS_DOT_H

/*
Authors: Sarah Al-Assam, Stephen Clark, Dieter Jaksch (Oxford) and Chris Goodyer (NAG)
$LastChangedDate: 2015-03-15 14:56:09 +0000 (Sun, 15 Mar 2015) $
(c) University of Oxford 2012-2015
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

/*! Reserved IDs for marking the start and end of networks */
#define TNT_NWSTART_ID 1
#define TNT_NWEND_ID 2

#define TNT_ID_OFFS 3 /* Number to add to all IDs to make sure that they do no interfere with the reserved ones */
#define TNT_MAXLEGNUM 80 /* Max id. Here the number is chosen as the difference between the min and max alphanumeric character codes with the offset added (+ a small buffer) */

/*! Fixed length for all strings. */
#define TNT_STRLEN 1000

/*! Maximum number of quantum numbers to print */
#define MAX_QN_PRINT 10

/*! Length for error messages */
#define TNT_MSG_LEN 300

/*! Define the width of messages to screen */
#define TNT_TEXT_WIDTH 524

/*! Maximum number of nodes that can be contracted in one function call */
#define TNT_MAXNUMNODES 12
#define TNT_MAXNUMLOAD 12

/*! Code for currently allowed symmetry type for the system - U(1). */
#define TNT_SYMM_NONE 0
#define TNT_SYMM_U1 1

/*! Code for currently allowed SVD routine types */
#define TNT_LAPACK_SDD 1
#define TNT_LAPACK_SVD 2
#define TNT_MAGMA_SDD 3

/*! Code for general preference between currently allowed LAPACK-like linear algebra libraries */
#define TNT_USE_LAPACK 0
#define TNT_USE_MAGMA 1

/*! Code for currently supported ARPACK libraries */
#define TNT_NAG_EIGS 1
#define TNT_ARPACK_EIGS 2

/*! Directions of legs for QN conservation */
#define TNT_QN_IN 1
#define TNT_QN_OUT -1

/*! Default tolerance for assuming values are zero. */
#define TNT_DEFAULT_TOL 1e-14

/*! Default maximum ratio between imaginary and real part of complex numbers to assume values are real. */
#define TNT_COMP_TOL 1e-8

/*! Define function for setting the random number seed */
#define tntSysRandSeedSet(A) srand48(A)
#define tntRandSeedSet(A) srand48((long) A)
#define tntRand() drand48()

/* Define the MIN and MAX macros. */
#if !defined(MAX)
#define MAX(A, B) ((A) > (B) ? (A): (B))
#endif

#if !defined(MIN)
#define MIN(A, B) ((A) < (B) ? (A): (B))
#endif

typedef const char * tntLegLabel;

/*! \ingroup node 
 * \brief Definition of type for network nodes. */
typedef struct tnode *tntNode;

/*! \ingroup publicnetwork 
 * \brief Definition of type for network. */
typedef struct tnetwork *tntNetwork;

/*! Array for integer values.*/
typedef struct int_arr {
    int *vals;  /*!< \brief Pointer to the first element in the integer array. */
    unsigned sz; /*!< \brief The number of elements in the array. */
    unsigned numrows; /*!< \brief The number of rows in the matrix. */
    unsigned numcols; /*!< \brief The number of columns in the matrix. */
} tntIntArray;

/*! Array for double values. */
typedef struct dbl_arr {
    double *vals; /*!< \brief Pointer to the first element in the double array */
    unsigned sz; /*!< \brief The number of elements in the array. */
    unsigned numrows; /*!< \brief The number of rows in the matrix. */
    unsigned numcols; /*!< \brief The number of columns in the matrix. */
} tntDoubleArray;

/*! Array for complex values. */
typedef struct cmp_arr {
    tntComplex *vals; /*!< \brief Pointer to the first element in the tntComplex array */
    unsigned sz; /*!< \brief The number of elements in the array. */
    unsigned numrows; /*!< \brief The number of rows in the matrix. */
    unsigned numcols; /*!< \brief The number of columns in the matrix. */
} tntComplexArray;

/*! Array for tntNodes. */
typedef struct tn_arr {
    tntNode *vals; /*!< \brief Pointer to the first element in the tntNode array */
    unsigned sz; /*!< \brief The number of elements in the array. */
} tntNodeArray;

/*! Array for strings (maximum length TNT_STRLEN). */
typedef struct str_arr {
    char **vals; /*!< \brief Pointer to the first element in the array of strings */
    unsigned sz; /*!< \brief The number of elements in the array. */
} tntStringArray;


/*! \ingroup core
 * Structure that holds all the operators required to calculate expectation values in 1D systems */
typedef struct expec_op {
    tntNodeArray ExOpOs;    /*!< \brief Array that contains the single site expectation values to calculate. */
    tntStringArray LbOpOs;  /*!< \brief Array that contains strings for labelling the above operators. */
    tntNodeArray ExOp2nn;   /*!< \brief Array that contains the nearest-neighbour expectation values to calculate. */
    tntStringArray LbOp2nn; /*!< \brief Array that contains strings for labelling the above operators. */
    tntNodeArray ExOp2cs;   /*!< \brief Array that contains the centre-to-site expectation values to calculate.  */
    tntStringArray LbOp2cs; /*!< \brief Array that contains strings for labelling the above operators. */
    tntNodeArray ExOp2ap;   /*!< \brief Array that contains the all pairs of two-site expectation values to calculate.*/
    tntStringArray LbOp2ap; /*!< \brief Array that contains strings for labelling the above operators. */
} tntExOp;

/*! Definition of type for defining the system. There is one global variable of this type used in the library functions. */
typedef struct tsystem tntSystem;

#endif
