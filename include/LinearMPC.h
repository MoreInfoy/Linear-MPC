//
// Created by nimapng on 6/7/21.
//

#ifndef LINEARMPC_LINEARMPC_H
#define LINEARMPC_LINEARMPC_H


#include <iostream>
#include <vector>
#include "qpOASES.hpp"
#include "MatrixDefinition.h"

using namespace std;
using namespace qpOASES;

class LinearMPC {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
public:
    LinearMPC();

    void reset();

    void setAkBk(Ref<Mat_Ak> Ak, Ref<Mat_Bk> Bk);

    void setAkBk(Ref<Mat_Ak> Ak, Ref<Mat_Bk> Bk, size_t k);

    void setStateConstraints(Ref<Mat_Ccx> Ccx, Ref<Vec_bcx> bcx);

    void setStateConstraints(Ref<Mat_Ccx> Ccx, Ref<Vec_bcx> bcx, size_t k);

    void setInputConstraints(Ref<Mat_Ccu> Ccu, Ref<Vec_bcu> bcu);

    void setInputConstraints(Ref<Mat_Ccu> Ccu, Ref<Vec_bcu> bcu, size_t k);

    void setInputBounds(Ref<Vec_U_Bounds> lb, Ref<Vec_U_Bounds> ub);

    void setInputBounds(Ref<Vec_U_Bounds> lb, Ref<Vec_U_Bounds> ub, size_t k);

    void setWeightMatrix(Ref<Mat_Qx> Qx, Ref<Mat_Ru> Ru);

    void setWeightMatrix(Ref<Mat_Qx> Qx, Ref<Mat_Ru> Ru, size_t k);

    void setInitialStateAndRef(Ref<Vec> x0, Ref<Vec> x_ref);

    ConstVecRef getSolution();

    Vec getOptimalTraj();

    real_t getCost();

private:
    QProblem solver;
    Mat_Ak Ak_vec[HORIZON];
    Mat_Bk Bk_vec[HORIZON];
    Mat_Ccx Ccx_vec[HORIZON];
    Mat_Ccu Ccu_vec[HORIZON];
    Vec_bcx bcx_vec[HORIZON];
    Vec_bcu bcu_vec[HORIZON];
    Vec_U_Bounds lb_vec[HORIZON];
    Vec_U_Bounds ub_vec[HORIZON];
    Mat_Qx Qx_vec[HORIZON];
    Mat_Ru Ru_vec[HORIZON];

    Vec g, c, bx, bu, lb, ub, sol, x0, x_ref;
    Mat_qp H, C, Sx, Su, Cx, Cu, Q, R;

    void solve();

    void computeSxSu();

    void computeCxbx();

    void computeCubu();

    void computelbub();

    void computeQR();
};


#endif //LINEARMPC_LINEARMPC_H
