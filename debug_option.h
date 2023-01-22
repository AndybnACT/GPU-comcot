#ifndef DBG_OPT
#define DBG_OPT

/* #define DEBUG_ALL_GRID */
/* #define DEBUG_CORE */

/* #define RELDIFF */

#ifdef RELDIFF
#define DIFFRATE 0.001
#else
#define TOLERANCE 1.0e-4
#endif /* RELDIFF */

#ifdef DEBUG_CORE
#define DEBUG_CORE_MASS
#define DEBUG_CORE_MOMT
#define DEBUG_CORE_JNQ
/* #define DEBUG_CORE_JNZ */
#endif /* DEBUG_CORE  */

#endif /* !DBG_OPT */
