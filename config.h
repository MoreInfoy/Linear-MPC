//
// Created by nimapng on 6/7/21.
//

#ifndef LINEAR_MPC_CONFIG_H
#define LINEAR_MPC_CONFIG_H

#define real_t double

#define Ak_ROWS 2
#define Bk_COLS 1

#define WITH_STATE_CONSTRAINTS

#ifdef WITH_STATE_CONSTRAINTS
#define SC_ROWS 1
#else
#define SC_ROWS 0
#endif

//#define WITH_OUTPUT
//
//#ifdef WITH_OUTPUT
//#define Ck_ROWS 12
//#define Dk_COLS 12
//
//#define WITH_OUTPUT_CONSTRAINTS
//#ifdef WITH_OUTPUT_CONSTRAINTS
//#define OC_ROWS 12
//#endif

//#else
//#define Ck_ROWS 0
//#define Dk_COLS 0
//#define OC_ROWS 0
//#endif

#define WITH_INPUT_BOUNDS

#define WITH_INPUT_CONSTRAINTS

#ifdef WITH_INPUT_CONSTRAINTS
#define UC_ROWS 1
#else
#define UC_ROWS 0
#endif

#define HORIZON 120

#endif //LINEAR_MPC_CONFIG_H
