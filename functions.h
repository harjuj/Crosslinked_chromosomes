#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "global.h"
#include "moves.h"

void normalize(std::vector<std::vector<double>> contacts) {
    double sum = 0;
    for (int i = 0; i < bin_num; i++) {
        for (int j = 0; j < bin_num; j++) {
            if (std::abs(i-j) > 1 && !(i==0 && j == bin_num - 1) && !(j==0 && i == bin_num - 1)) {
                sum += contacts[i][j];
            }
        }
    }
    for (int i = 0; i < bin_num; i++) {
        for (int j = 0; j < bin_num; j++) {
            contacts[i][j] *= float(bin_num) / sum;
        }
    }
}


#endif //FUNCTIONS_H
