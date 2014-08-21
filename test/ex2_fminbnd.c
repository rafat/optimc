/*
 ============================================================================
 Name        : optim.c
 Author      : Rafat
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include "../header/optimc.h"
#include "testfunctions.h"


int main(void) {
        double a,b,oup;

        a = 0.3;
        b = 1;

        custom_funcuni humps_min = {humps,0};

        oup = fminbnd(&humps_min,a,b);
        printf("OUP %g \n",oup);

        return 0;
}
