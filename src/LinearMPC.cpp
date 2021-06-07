//
// Created by nimapng on 6/7/21.
//

#include "LinearMPC.h"

LinearMPC::LinearMPC() : solver(Bk_COLS * HORIZON,
                                SC_ROWS + UC_ROWS, HST_SEMIDEF) {
    reset();
}

void LinearMPC::reset() {
    for (size_t i = 0; i < HORIZON; i++) {
        Ak_vec[i].setZero();
        Bk_vec[i].setZero();
        Ccx_vec[i].setZero();
        Ccu_vec[i].setZero();
        bcx_vec[i].setZero();
        bcu_vec[i].setZero();
        lb_vec[i].setZero();
        ub_vec[i].setZero();
        Qx_vec[i].setZero();
        Ru_vec[i].setZero();
        H.resize(Bk_COLS * HORIZON, Bk_COLS * HORIZON);
        C.resize(SC_ROWS + UC_ROWS, Bk_COLS * HORIZON);
        c.resize(SC_ROWS + UC_ROWS);
        lb.resize(Bk_COLS * HORIZON);
        ub.resize(Bk_COLS * HORIZON);
        g.resize(Bk_COLS * HORIZON);
        sol.resize(Bk_COLS * HORIZON);
    }

    Options opt;
    opt.setToMPC();
    opt.enableEqualities = BT_TRUE;
    solver.setOptions(opt);
}

void LinearMPC::setAkBk(Ref<Mat_Ak> Ak, Ref<Mat_Bk> Bk) {
    for (size_t i = 0; i < HORIZON; i++) {
        Ak_vec[i] = Ak;
        Bk_vec[i] = Ak;
    }
}


void LinearMPC::setAkBk(Ref<Mat_Ak> Ak, Ref<Mat_Bk> Bk, size_t k) {
    if (k > HORIZON) {
        throw std::runtime_error("k > HORIZON");
    }

    Ak_vec[k] = Ak;
    Bk_vec[k] = Ak;
}

void LinearMPC::setStateConstraints(Ref<Mat_Ccx> Ccx, Ref<Vec_bcx> bcx) {
    for (size_t k = 0; k < HORIZON; k++) {
        Ccx_vec[k] = Ccx;
        bcx_vec[k] = bcx;
    }
}

void LinearMPC::setStateConstraints(Ref<Mat_Ccx> Ccx, Ref<Vec_bcx> bcx, size_t k) {
    if (k > HORIZON) {
        throw std::runtime_error("k > HORIZON");
    }
    Ccx_vec[k] = Ccx;
    bcx_vec[k] = bcx;
}

void LinearMPC::setInputConstraints(Ref<Mat_Ccu> Ccu, Ref<Vec_bcu> bcu) {
    for (size_t k = 0; k < HORIZON; k++) {
        Ccu_vec[k] = Ccu;
        bcu_vec[k] = bcu;
    }
}

void LinearMPC::setInputConstraints(Ref<Mat_Ccu> Ccu, Ref<Vec_bcu> bcu, size_t k) {
    if (k > HORIZON) {
        throw std::runtime_error("k > HORIZON");
    }
    Ccu_vec[k] = Ccu;
    bcu_vec[k] = bcu;
}

void LinearMPC::setInputBounds(Ref<Vec_U_Bounds> lb, Ref<Vec_U_Bounds> ub) {
    for (size_t k = 0; k < HORIZON; k++) {
        lb_vec[k] = lb;
        ub_vec[k] = ub;
    }
}

void LinearMPC::setInputBounds(Ref<Vec_U_Bounds> lb, Ref<Vec_U_Bounds> ub, size_t k) {
    if (k > HORIZON) {
        throw std::runtime_error("k > HORIZON");
    }
    lb_vec[k] = lb;
    ub_vec[k] = ub;
}

void LinearMPC::setWeightMatrix(Ref<Mat_Qx> Qx, Ref<Mat_Ru> Ru) {
    for (size_t i = 0; i < HORIZON; i++) {
        Qx_vec[i] = Qx;
        Ru_vec[i] = Ru;
    }
}

void LinearMPC::setWeightMatrix(Ref<Mat_Qx> Qx, Ref<Mat_Ru> Ru, size_t k) {
    if (k >= HORIZON) {
        throw std::runtime_error("k > HORIZON");
    }
    Qx_vec[k] = Qx;
    Ru_vec[k] = Ru;
}


void LinearMPC::setInitialStateAndRef(Ref<Vec> x0, Ref<Vec> x_ref) {
    this->x0 = x0;
    this->x_ref = x_ref;
}

ConstVecRef LinearMPC::getSolution() {
    solve();
    return ConstVecRef(sol);
}

void LinearMPC::solve() {
    /* add formulation process */

    int_t nWSR = 1000;
    solver.init(H.data(), g.data(), C.data(), lb.data(), ub.data(), c.data(), NULL, nWSR);
    if(solver.getPrimalSolution(sol.data()) != SUCCESSFUL_RETURN)
    {
        throw std::runtime_error("qp solver failed");
    }

}





