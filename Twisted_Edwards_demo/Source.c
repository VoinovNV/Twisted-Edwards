#include <stdio.h>
#include "tests.h"
#include <assert.h>

int main() {
const int size=12;
void (*f[size])()={
    test_add,
    test_double,
    test_pow_montgomery,
    test_form_changing,
    test_triple,
    test_quintuple,
    test_pow_bin,
    test_pow_NAF,
    test_pow_NAF_2_w,
    test_pow_NAF_l_w,
    test_pow_NAF_mbw,
    ak_signkey_test_edwards
};
for (int i=0;i<size;i++) {
    printf("Test %d: \n",i);
    f[i]();
    printf("OK!\n");
}
printf("\nTIME:!\n");
    get_time();

    printf("\nEND!\n");
}
