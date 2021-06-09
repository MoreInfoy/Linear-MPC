//
// Created by nimapng on 6/7/21.
//

#ifndef LINEARMPC_MATRIXDEFINITION_H
#define LINEARMPC_MATRIXDEFINITION_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include "config.h"

using namespace Eigen;

namespace MPC {
    using Vec = typename Eigen::Matrix<real_t, Dynamic, 1>;

    using Mat = typename Eigen::Matrix<real_t, Dynamic, Dynamic>;

    using Mat_qp = typename Eigen::Matrix<real_t, Dynamic, Dynamic, RowMajor>;

    using Mat_Ak = typename Eigen::Matrix<real_t, Ak_ROWS, Ak_ROWS>;
    using Mat_Bk = typename Eigen::Matrix<real_t, Ak_ROWS, Bk_COLS>;

    using Mat_Qx = typename Eigen::Matrix<real_t, Ak_ROWS, Ak_ROWS>;
    using Mat_Ru = typename Eigen::Matrix<real_t, Bk_COLS, Bk_COLS>;

    using Mat_Ccx = typename Eigen::Matrix<real_t, SC_ROWS, Ak_ROWS>;
    using Vec_bcx = typename Eigen::Matrix<real_t, SC_ROWS, 1>;

    using Mat_Ccu = typename Eigen::Matrix<real_t, UC_ROWS, Bk_COLS>;
    using Vec_bcu = typename Eigen::Matrix<real_t, UC_ROWS, 1>;

#ifdef WITH_INPUT_BOUNDS
    using Vec_U_Bounds = typename Eigen::Matrix<real_t, Bk_COLS, 1>;
#else
    using Vec_U_Bounds = typename Eigen::Matrix<real_t, 0, 0>;
#endif

#ifdef WITH_OUTPUT
    using Mat_Ck = typename Eigen::Matrix<real_t, Ck_ROWS, Ak_ROWS, RowMajor>;
    using Mat_Dk = typename Eigen::Matrix<real_t, Ck_ROWS, Dk_COLS, RowMajor>;
#ifdef WITH_OUTPUT_CONSTRAINTS
    using Mat_Ccy = typename Eigen::Matrix<real_t, OC_ROWS, Ak_ROWS, RowMajor>;
    using Vec_bcy = typename Eigen::Matrix<real_t, OC_ROWS, 1, RowMajor>;
#endif
#endif

    typedef const Eigen::Ref<const Vec> ConstVecRef;

    typedef const Eigen::Ref<const Mat> ConstMatRef;

    typedef Eigen::Ref<Vec> VecRef;

    typedef Eigen::Ref<Mat> MatRef;
}

#endif //LINEARMPC_MATRIXDEFINITION_H
