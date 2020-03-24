#define LIBAKRYPT_HAVE_SYSTYPES_H
#include <stdio.h>
#include <time.h>
#include "tests.h"
#include <assert.h>
#include <ak_sign.h>

int main() {

	printf("Test 1: ");
	test_add();
	printf("done\nTest 2: ");
	test_double();
	printf("done\nTest 3: ");
	test_pow();
	printf("done\nTest 4: ");
	test_form_changing();
	printf("done\nTest 5: ");
	assert(ak_signkey_test_edwards());
	printf("done\n");
	int start, end;
	
	for (int j = 1; j <= 150; j+=50) {
		printf("j=%d\n", j);
		start = clock();
		for (int i = 0; i < j; i++) {
			ak_signkey_test_edwards();
		}
		end = clock();
		printf("Edwards: %lf\n", (float)(end - start) / (CLOCKS_PER_SEC));
		start = clock();
		for (int i = 0; i < j; i++) {
			my_ak_signkey_test();
		}
		end = clock();
		printf("Weierstrass: %lf\n", (float)(end - start) / (CLOCKS_PER_SEC));
		
	}
    
    test_triple();

    test_quintuple();
	printf("OK!\n");

}
