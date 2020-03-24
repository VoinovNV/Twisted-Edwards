#pragma once
#include "ecurve.h"
#include <ak_sign.h>

void ak_signkey_context_sign_const_values_edwards(ak_signkey sctx, 
	ak_uint64 *k, ak_uint64 *e, ak_pointer out, ak_ecurve ec);

bool_t ak_verifykey_context_verify_hash_edwards(ak_verifykey pctx,
	const ak_pointer hash, const size_t hsize, ak_pointer sign, ak_ecurve ec);
