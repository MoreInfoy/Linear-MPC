//
// Created by nimapng on 6/7/21.
//

#include "LinearMPC.h"

int main()
{
    LinearMPC linearMpc;
    MPC::Mat_Ak Ak;
    MPC::Mat_Bk Bk;
    Ak << 1.0, 0.05, 0, 1.0;
    Bk << 0.0012, 0.05;
    linearMpc.setAkBk(Ak, Bk);

    MPC::Mat_Ccx Ccx;
    Ccx << 1.0, 0;
    MPC::Vec_bcx cxlb,cxub;
    cxlb << -3.0;
    cxub << 3.0;
    linearMpc.setStateConstraints(Ccx, cxlb, cxub);

    MPC::Mat_Ccu Ccu;
    Ccu << 0.5;
    MPC::Vec_bcu culb,cuub;
    culb << -3.0;
    cuub << 3.0;
    linearMpc.setInputConstraints(Ccu, culb, cuub);

    MPC::Vec_U_Bounds lb, ub;
    lb << -5.5;
    ub << 4.6;
    linearMpc.setInputBounds(lb, ub);

    MPC::Mat_Qx Q;
    Q << 10.0, 0, 0, 10.0;
    MPC::Mat_Ru R;
    R << 0.1;
    linearMpc.setWeightMatrix(Q, R);

    MPC::Vec x0, x_ref;
    x0.resize(2);
    x0 << 2.0, 0.0;
    x_ref.resize(HORIZON * Ak_ROWS);
    x_ref.setZero();

    linearMpc.setInitialStateAndRef(x0, x_ref);

    std::cout << "sol: " << linearMpc.getSolution().transpose() << std::endl;

    linearMpc.outputAllDataToFile("data.txt");

    return 0;
}