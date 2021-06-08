//
// Created by nimapng on 6/7/21.
//

#include "LinearMPC.h"

int main()
{
    LinearMPC linearMpc;
    linearMpc.reset();
    linearMpc.outputAllDataToFile("data.txt");

    return 0;
}