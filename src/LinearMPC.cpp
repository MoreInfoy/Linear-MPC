//
// Created by nimapng on 6/7/21.
//

#include "LinearMPC.h"

LinearMPC::LinearMPC() : solver(Bk_COLS * HORIZON,
                                (SC_ROWS + UC_ROWS) * HORIZON, HST_SEMIDEF) {
    reset();
}

void LinearMPC::reset() {
    for (size_t i = 0; i < HORIZON; i++) {
        Ak_vec[i].setZero();
        Bk_vec[i].setZero();
        Ccx_vec[i].setZero();
        Ccu_vec[i].setZero();
        lbcx_vec[i].setZero();
        ubcx_vec[i].setZero();
        lbcu_vec[i].setZero();
        ubcu_vec[i].setZero();
        lb_vec[i].setZero();
        ub_vec[i].setZero();
        Qx_vec[i].setZero();
        Ru_vec[i].setZero();
    }
    H.resize(Bk_COLS * HORIZON, Bk_COLS * HORIZON);
    H.setZero();
    C.resize((SC_ROWS + UC_ROWS) * HORIZON, Bk_COLS * HORIZON);
    C.setZero();
    Sx.resize(Ak_ROWS * HORIZON, Ak_ROWS);
    Sx.setZero();
    Su.resize(Ak_ROWS * HORIZON, Bk_COLS * HORIZON);
    Su.setZero();
    Cx.resize(SC_ROWS * HORIZON, Ak_ROWS * HORIZON);
    Cx.setZero();
    lbx.resize(SC_ROWS * HORIZON);
    lbx.setZero();
    ubx.resize(SC_ROWS * HORIZON);
    ubx.setZero();
    Cu.resize(UC_ROWS * HORIZON, Bk_COLS * HORIZON);
    Cu.setZero();
    lbu.resize(UC_ROWS * HORIZON);
    lbu.setZero();
    ubu.resize(UC_ROWS * HORIZON);
    ubu.setZero();
    clb.resize((SC_ROWS + UC_ROWS) * HORIZON);
    clb.setZero();
    cub.resize((SC_ROWS + UC_ROWS) * HORIZON);
    cub.setZero();
    lb.resize(Bk_COLS * HORIZON);
    lb.setZero();
    ub.resize(Bk_COLS * HORIZON);
    ub.setZero();
    Q.resize(Ak_ROWS * HORIZON, Ak_ROWS * HORIZON);
    Q.setZero();
    R.resize(Bk_COLS * HORIZON, Bk_COLS * HORIZON);
    R.setZero();
    g.resize(Bk_COLS * HORIZON);
    g.setZero();
    sol.resize(Bk_COLS * HORIZON);
    sol.setZero();

    x_ref.resize(Ak_ROWS * HORIZON);
    x_ref.setZero();
    x0.resize(Ak_ROWS);
    x0.setZero();

    Options opt;
    opt.setToMPC();
    opt.enableEqualities = BT_TRUE;
    solver.setOptions(opt);
}

void LinearMPC::setAkBk(Ref<Mat_Ak> Ak, Ref<Mat_Bk> Bk) {
    for (size_t i = 0; i < HORIZON; i++) {
        Ak_vec[i] = Ak;
        Bk_vec[i] = Bk;
    }
}


void LinearMPC::setAkBk(Ref<Mat_Ak> Ak, Ref<Mat_Bk> Bk, size_t k) {
    if (k > HORIZON) {
        throw std::runtime_error("k > HORIZON");
    }

    Ak_vec[k] = Ak;
    Bk_vec[k] = Bk;
}

void LinearMPC::setStateConstraints(Ref<Mat_Ccx> Ccx, Ref<Vec_bcx> lbcx, Ref<Vec_bcx> ubcx) {
    for (size_t k = 0; k < HORIZON; k++) {
        Ccx_vec[k] = Ccx;
        lbcx_vec[k] = lbcx;
        ubcx_vec[k] = ubcx;
    }
}

void LinearMPC::setStateConstraints(Ref<Mat_Ccx> Ccx, Ref<Vec_bcx> lbcx, Ref<Vec_bcx> ubcx, size_t k) {
    if (k > HORIZON) {
        throw std::runtime_error("k > HORIZON");
    }
    Ccx_vec[k] = Ccx;
    lbcx_vec[k] = lbcx;
    ubcx_vec[k] = ubcx;
}

void LinearMPC::setInputConstraints(Ref<Mat_Ccu> Ccu, Ref<Vec_bcu> lbcu, Ref<Vec_bcu> ubcu) {
    for (size_t k = 0; k < HORIZON; k++) {
        Ccu_vec[k] = Ccu;
        lbcu_vec[k] = lbcu;
        ubcu_vec[k] = ubcu;
    }
}

void LinearMPC::setInputConstraints(Ref<Mat_Ccu> Ccu, Ref<Vec_bcu> lbcu, Ref<Vec_bcu> ubcu, size_t k) {
    if (k > HORIZON) {
        throw std::runtime_error("k > HORIZON");
    }
    Ccu_vec[k] = Ccu;
    lbcu_vec[k] = lbcu;
    ubcu_vec[k] = ubcu;
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
    computeSxSu();
    computeCxbx();
    computeCubu();
    computelbub();
    computeQR();

    C.topRows(SC_ROWS * HORIZON).noalias() = Cx * Su;
    C.bottomRows(UC_ROWS * HORIZON) = Cu;

    clb.head(SC_ROWS * HORIZON) = lbx - Cx * Sx * x0;
    clb.tail(UC_ROWS * HORIZON) = lbu;

    cub.head(SC_ROWS * HORIZON) = ubx - Cx * Sx * x0;
    cub.tail(UC_ROWS * HORIZON) = ubu;

    H = R + Su.transpose() * Q * Su;
    g = Su.transpose() * Q * Sx * x0 - Su.transpose() * Q * x_ref;

    int_t nWSR = 1000;
    solver.reset();
    Options opt;
    opt.setToMPC();
    opt.enableEqualities = BT_TRUE;
    solver.setOptions(opt);
    solver.init(H.data(), g.data(), C.data(), lb.data(), ub.data(), clb.data(), cub.data(), nWSR);
    if (solver.isSolved()) {
        solver.getPrimalSolution(sol.data());
    } else {
        throw std::runtime_error("qp solver failed");
    }

}

void LinearMPC::computeSxSu() {
    for (int k = 0; k < HORIZON; k++) {
        Sx.middleRows(k * Ak_ROWS, Ak_ROWS) = Ak_vec[k];
        if (k == 0) {
            Su.topLeftCorner<Ak_ROWS, Bk_COLS>() = Bk_vec[k];
        } else {
            Su.block(k * Ak_ROWS, 0, Ak_ROWS, k * Bk_COLS)
                    = Ak_vec[k] * (Su.block((k - 1) * Ak_ROWS, 0, Ak_ROWS, k * Bk_COLS));
            Su.block<Ak_ROWS, Bk_COLS>(k * Ak_ROWS, k * Bk_COLS) = Bk_vec[k];
        }
    }
}

void LinearMPC::computeCxbx() {
    for (int k = 0; k < HORIZON; k++) {
        Cx.block<SC_ROWS, Ak_ROWS>(k * SC_ROWS, k * Ak_ROWS) = Ccx_vec[k];
        lbx.block<SC_ROWS, 1>(k * SC_ROWS, 0) = lbcx_vec[k];
        ubx.block<SC_ROWS, 1>(k * SC_ROWS, 0) = ubcx_vec[k];
    }
}

void LinearMPC::computeCubu() {
    for (int k = 0; k < HORIZON; k++) {
        Cu.block<UC_ROWS, Bk_COLS>(k * UC_ROWS, k * Bk_COLS) = Ccu_vec[k];
        lbu.block<UC_ROWS, 1>(k * UC_ROWS, 0) = lbcu_vec[k];
        ubu.block<UC_ROWS, 1>(k * UC_ROWS, 0) = ubcu_vec[k];
    }
}

void LinearMPC::computelbub() {
    for (int k = 0; k < HORIZON; k++) {
        lb.block<Bk_COLS, 1>(k * Bk_COLS, 0) = lb_vec[k];
        ub.block<Bk_COLS, 1>(k * Bk_COLS, 0) = ub_vec[k];
    }
}

void LinearMPC::computeQR() {
    for (int k = 0; k < HORIZON; k++) {
        Q.block<Ak_ROWS, Ak_ROWS>(k * Ak_ROWS, k * Ak_ROWS) = Qx_vec[k];
        R.block<Bk_COLS, Bk_COLS>(k * Bk_COLS, k * Bk_COLS) = Ru_vec[k];
    }
}

Vec LinearMPC::getOptimalTraj() {
    return Sx * x0 + Su * sol;
}

real_t LinearMPC::getCost() {
    if (solver.isSolved()) {
        return solver.getObjVal();
    } else {
        throw std::runtime_error("qp solver failed, cost can not be obtained");
    }
}

void LinearMPC::outputAllDataToFile(string file_name) {
    ofstream outfile(file_name);
    if (!outfile.is_open()) {
        throw std::runtime_error("[LinearMPC::outputAllDataToFile] The file can not be opened");
    }
    outfile << "-----------------------H------------------------" << endl
            << H << endl
            << "-----------------------g------------------------" << endl
            << g.transpose() << endl
            << "-----------------------C------------------------" << endl
            << C << endl
            << "-----------------------c_lb------------------------" << endl
            << clb.transpose() << endl
            << "-----------------------c_ub------------------------" << endl
            << cub.transpose() << endl
            << "-----------------------lb------------------------" << endl
            << lb.transpose() << endl
            << "-----------------------ub------------------------" << endl
            << ub.transpose() << endl
            << "-----------------------sol------------------------" << endl
            << sol.transpose() << endl
            << "-----------------------optimal trajectory----------------------" << endl
            << (Sx * x0 + Su * sol).transpose() << endl;
    outfile.close();
}






