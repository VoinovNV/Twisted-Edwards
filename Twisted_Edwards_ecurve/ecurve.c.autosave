#include "ecurve.h"
#include <string.h>  

void ak_epoint_set_epoint(ak_epoint ep1, ak_epoint ep2, ak_ecurve ec)
{
	memcpy(ep1->X, ep2->X, ec->size * sizeof(ak_uint64));
	memcpy(ep1->Y, ep2->Y, ec->size * sizeof(ak_uint64));
	memcpy(ep1->Z, ep2->Z, ec->size * sizeof(ak_uint64));
}

void ak_epoint_set(ak_epoint ep, ak_ecurve ec) {
	memcpy(ep->X, ec->point.X, ec->size * sizeof(ak_uint64));
	memcpy(ep->Y, ec->point.Y, ec->size * sizeof(ak_uint64));
	memcpy(ep->Z, ec->point.Z, ec->size * sizeof(ak_uint64));
}

void ak_epoint_set_as_unit(ak_epoint ep, ak_ecurve ec)
{
	ak_mpzn_set_ui(ep->X, ec->size, 0);
	ak_mpzn_set_ui(ep->Y, ec->size, 1);
	ak_mpzn_set_ui(ep->Z, ec->size, 1);
}

void ak_epoint_reduce(ak_epoint ep, ak_ecurve ec)
{
	if (ak_mpzn_cmp_ui(ep->Z, ec->size, 0)) return;
	ak_mpznmax u, one = ak_mpznmax_one;
	ak_mpzn_set_ui(u, ec->size, 2);
	ak_mpzn_sub(u, ec->p, u, ec->size);
	ak_mpzn_modpow_montgomery(u, ep->Z, u, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(u, u, one, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(ep->X, ep->X, u, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(ep->Y, ep->Y, u, ec->p, ec->n, ec->size);
	ak_mpzn_set_ui(ep->Z, ec->size, 1);
}

bool_t ak_epoint_is_ok(ak_epoint ep, ak_ecurve ec) {
	ak_mpznmax a, b, c, d;
	ak_mpzn_mul_montgomery(a, ep->X, ep->X, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(d, a, ec->e, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(b, ep->Y, ep->Y, ec->p, ec->n, ec->size);
	ak_mpzn_add_montgomery(d, d, b, ec->p, ec->size);
	ak_mpzn_mul_montgomery(c, ep->Z, ep->Z, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(d, d, c, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(c, c, c, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(a, a, b, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(a, a, ec->d, ec->p, ec->n, ec->size);
	ak_mpzn_add_montgomery(a, a, c, ec->p, ec->size);
	if (ak_mpzn_cmp(a, d, ec->size)) return ak_false;
	return ak_true;
}
void ak_epoint_double(ak_epoint ep, ak_ecurve ec) {
	ak_mpznmax b, c, d, e, f, h, j, st;
	ak_mpzn_add_montgomery(b, ep->X, ep->Y, ec->p, ec->size);
	ak_mpzn_mul_montgomery(b, b, b, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(c, ep->X, ep->X, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(d, ep->Y, ep->Y, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(e, ec->e, c, ec->p, ec->n, ec->size);
	ak_mpzn_add_montgomery(f, e, d, ec->p, ec->size);
	ak_mpzn_mul_montgomery(h, ep->Z, ep->Z, ec->p, ec->n, ec->size);
	ak_mpzn_lshift_montgomery(j, h, ec->p, ec->size);
	ak_mpzn_sub(st, ec->p, j, ec->size);
	ak_mpzn_add_montgomery(j, f, st, ec->p, ec->size);
	ak_mpzn_sub(st, ec->p, c, ec->size);
	ak_mpzn_add_montgomery(b, b, st, ec->p, ec->size);
	ak_mpzn_sub(st, ec->p, d, ec->size);
	ak_mpzn_add_montgomery(b, b, st, ec->p, ec->size);
	ak_mpzn_mul_montgomery(ep->X, b, j, ec->p, ec->n, ec->size);
	ak_mpzn_sub(st, ec->p, d, ec->size);
	ak_mpzn_add_montgomery(e, e, st, ec->p, ec->size);
	ak_mpzn_mul_montgomery(ep->Y, f, e, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(ep->Z, f, j, ec->p, ec->n, ec->size);
}
void ak_epoint_add(ak_epoint ep1, ak_epoint ep2, ak_ecurve ec) {
	ak_mpznmax a, b, c, d, e, f, g, st;
	ak_mpzn_mul_montgomery(a, ep1->Z, ep2->Z, ec->p, ec->n, ec->size); 
	ak_mpzn_mul_montgomery(b, a, a, ec->p, ec->n, ec->size); 
	ak_mpzn_mul_montgomery(c, ep1->X, ep2->X, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(d, ep1->Y, ep2->Y, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(e, c, d, ec->p, ec->n, ec->size); 
	ak_mpzn_mul_montgomery(e, e, ec->d, ec->p, ec->n, ec->size); 
	ak_mpzn_sub(st, ec->p, e, ec->size);
	ak_mpzn_add_montgomery(f, b, st, ec->p, ec->size);
	ak_mpzn_add_montgomery(g, b, e, ec->p, ec->size); 
	ak_mpzn_add_montgomery(e, ep1->X, ep1->Y, ec->p, ec->size);
	ak_mpzn_add_montgomery(b, ep2->X, ep2->Y, ec->p, ec->size);
	ak_mpzn_mul_montgomery(b, b, e, ec->p, ec->n, ec->size); 
	ak_mpzn_sub(st, ec->p, c, ec->size);
	ak_mpzn_add_montgomery(b, b, st, ec->p, ec->size);
	ak_mpzn_sub(st, ec->p, d, ec->size);
	ak_mpzn_add_montgomery(b, b, st, ec->p, ec->size);
	ak_mpzn_mul_montgomery(b, b, f, ec->p, ec->n, ec->size); 
	ak_mpzn_mul_montgomery(ep1->X, b, a, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(b, ec->e, c, ec->p, ec->n, ec->size);
	ak_mpzn_sub(st, ec->p, b, ec->size);
	ak_mpzn_add_montgomery(b, d, st, ec->p, ec->size);
	ak_mpzn_mul_montgomery(b, b, g, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(ep1->Y, b, a, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(ep1->Z, f, g, ec->p, ec->n, ec->size);
}


void ak_epoint_pow(ak_epoint ep1, ak_epoint ep2, ak_uint64 * a, size_t size, ak_ecurve ec) {
	ak_uint64 uk = 0;
	long long int i, j;
	struct epoint Q, R; 
	ak_epoint_set_as_unit(&Q, ec);
	ak_epoint_set_epoint(&R, ep2, ec);
	for (i = size - 1; i >= 0; i--) {
		uk = a[i];
		for (j = 0; j < 64; j++) {
			if (uk & 0x8000000000000000LL) {ak_epoint_add(&Q, &R, ec); ak_epoint_double(&R, ec);}
			else {	ak_epoint_add(&R, &Q, ec); ak_epoint_double(&Q, ec);}
			uk <<= 1;
		}
	}
	ak_epoint_set_epoint(ep1, &Q, ec);
}


void ak_epoint_to_wpoint(ak_epoint ep, ak_ecurve ec, ak_wpoint wp, ak_wcurve wc) {
	if (ak_mpzn_cmp_ui(ep->Z, ec->size, 0)) {
		ak_wpoint_set_as_unit(wp, wc);
		return;
	}
	ak_mpznmax a_;
	ak_mpzn_sub(wp->z, ec->p, ep->Y, ec->size);
	ak_mpzn_add_montgomery(wp->z, wp->z, ep->Z, ec->p, ec->size);
	ak_mpzn_mul_montgomery(wp->z, wp->z, ep->X, ec->p, ec->n, ec->size);
	ak_mpzn_add_montgomery(a_, ep->Y, ep->Z, ec->p, ec->size);
	ak_mpzn_mul_montgomery(a_, a_, ec->S, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(wp->y, a_, ep->Z, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(wp->x, a_, ep->X, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(a_, ec->T, wp->z, ec->p, ec->n, ec->size);
	ak_mpzn_add_montgomery(wp->x,wp->x, a_, ec->p, ec->size);
}
void ak_wpoint_to_epoint(ak_wpoint wp, ak_wcurve wc, ak_epoint ep, ak_ecurve ec) {
	if (ak_mpzn_cmp_ui(wp->y, wc->size, 0)) {
		ak_epoint_set_as_unit(ep, ec);
		return;
	}
	ak_mpznmax a_,b_;
	ak_mpzn_mul_montgomery(a_, wp->z, ec->T, ec->p, ec->n, ec->size);
	ak_mpzn_sub(a_, wc->p, a_, ec->size);
	ak_mpzn_add_montgomery(a_, wp->x, a_, ec->p, ec->size);
	ak_mpzn_mul_montgomery(b_, wp->z, ec->S, ec->p, ec->n, ec->size);
	ak_mpzn_add_montgomery(ep->X, b_, a_, ec->p, ec->size);
	ak_mpzn_mul_montgomery(ep->Z, ep->X, wp->y, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(ep->X, a_, ep->X, ec->p, ec->n, ec->size);
	ak_mpzn_sub(ep->Y, ec->p, b_, ec->size);
	ak_mpzn_add_montgomery(ep->Y, ep->Y, a_, ec->p, ec->size);
	ak_mpzn_mul_montgomery(ep->Y, ep->Y, wp->y, ec->p, ec->n, ec->size);
}
void ak_epoint_triple(ak_epoint ep, ak_ecurve ec){
    ak_mpznmax a,b,c,d,e,f,g;
    ak_mpzn_mul_montgomery(a, ep->Y, ep->Y, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(b, ep->X, ep->X, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(b, b, ec->e, ec->p, ec->n, ec->size);
    ak_mpzn_add_montgomery(c, a, b, ec->p, ec->size);
    ak_mpzn_mul_montgomery(d, ep->Z, ep->Z, ec->p, ec->n, ec->size);
    ak_mpzn_lshift_montgomery(d,d,ec->p, ec->size);
    ak_mpzn_sub(e,ec->p,c,ec->size);
    ak_mpzn_add_montgomery(d, d, e, ec->p, ec->size);
    ak_mpzn_lshift_montgomery(d,d,ec->p, ec->size); //b=2(2*Z^2...)
    ak_mpzn_mul_montgomery(f, b, d, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(g, a, d, ec->p, ec->n, ec->size);
    ak_mpzn_sub(e,ec->p,b,ec->size);
    ak_mpzn_add_montgomery(e, a, e, ec->p, ec->size);
    ak_mpzn_mul_montgomery(e, e, c, ec->p, ec->n, ec->size);
    ak_mpzn_sub(a,ec->p,g,ec->size);
    ak_mpzn_add_montgomery(a, e, a, ec->p, ec->size);
    ak_mpzn_add_montgomery(b, e, f, ec->p, ec->size);
    ak_mpzn_add_montgomery(c, g, e, ec->p, ec->size);
    ak_mpzn_mul_montgomery(ep->X, ep->X, c, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(ep->X, ep->X, a, ec->p, ec->n, ec->size);
    ak_mpzn_sub(c,ec->p,e,ec->size);
    ak_mpzn_add_montgomery(c, c, f, ec->p, ec->size);
    ak_mpzn_mul_montgomery(ep->Y, ep->Y, c, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(ep->Y, ep->Y, b, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(ep->Z, ep->Z, a, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(ep->Z, ep->Z, b, ec->p, ec->n, ec->size);
}
void ak_epoint_quintuple(ak_epoint ep, ak_ecurve ec){
    ak_mpznmax x5,y5,z5,t1,t2,t5,t_;
    /*
    ak_mpzn_mul_montgomery(a, ep->Y, ep->Y, ec->p, ec->n, ec->size);
    ak_mpzn_lshift_montgomery(d,d,ec->p, ec->size);
    ak_mpzn_sub(e,ec->p,c,ec->size);
    ak_mpzn_add_montgomery(d, d, e, ec->p, ec->size);
    */
    ak_mpzn_mul_montgomery(z5, ep->Y, ep->Y, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(y5, ep->X, ep->X, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(y5, y5, ec->e, ec->p, ec->n, ec->size);
    ak_mpzn_add_montgomery(x5, z5, y5, ec->p, ec->size);
    ak_mpzn_mul_montgomery(t1, ep->Z, ep->Z, ec->p, ec->n, ec->size);
    ak_mpzn_lshift_montgomery(t1,t1,ec->p, ec->size);
    ak_mpzn_sub(t1,ec->p,t1,ec->size);
    ak_mpzn_add_montgomery(t1, t1, x5, ec->p, ec->size);
    ak_mpzn_lshift_montgomery(z5,z5,ec->p, ec->size);
    ak_mpzn_lshift_montgomery(y5,y5,ec->p, ec->size);
    ak_mpzn_mul_montgomery(z5, z5, t1, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(t1, t1, y5, ec->p, ec->n, ec->size);
    ak_mpzn_sub(y5,ec->p,y5,ec->size);
    ak_mpzn_add_montgomery(y5, y5, x5, ec->p, ec->size);
    ak_mpzn_mul_montgomery(y5, x5, y5, ec->p, ec->n, ec->size);
    ak_mpzn_add_montgomery(x5, y5, z5, ec->p, ec->size);
    ak_mpzn_sub(t2,ec->p,z5,ec->size);
    ak_mpzn_add_montgomery(t2, y5, t2, ec->p, ec->size);
    ak_mpzn_mul_montgomery(x5, x5, t2, ec->p, ec->n, ec->size);
    ak_mpzn_sub(t2,ec->p,t1,ec->size);
    ak_mpzn_add_montgomery(t2, y5, t2, ec->p, ec->size);
    ak_mpzn_add_montgomery(y5, y5, t1, ec->p, ec->size);
    ak_mpzn_mul_montgomery(y5, y5, t2, ec->p, ec->n, ec->size);
    ak_mpzn_add_montgomery(t2, t2, t1, ec->p, ec->size);
    ak_mpzn_mul_montgomery(t1, t1, x5, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(z5, z5, y5, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(x5, x5, t2, ec->p, ec->n, ec->size);
    ak_mpzn_sub(x5,ec->p,x5,ec->size);
    ak_mpzn_mul_montgomery(y5, t2, y5, ec->p, ec->n, ec->size);
    ak_mpzn_add_montgomery(t2, x5, z5, ec->p, ec->size);
    ak_mpzn_sub(t_,ec->p,z5,ec->size);
    ak_mpzn_add_montgomery(x5, x5, t_, ec->p, ec->size);
    ak_mpzn_mul_montgomery(x5, t2, x5, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(ep->X, ep->X, x5, ec->p, ec->n, ec->size);
    ak_mpzn_add_montgomery(z5, y5, t1, ec->p, ec->size);
    ak_mpzn_sub(t_,ec->p,t1,ec->size);
    ak_mpzn_add_montgomery(y5, y5, t_, ec->p, ec->size);
    ak_mpzn_mul_montgomery(y5, y5, z5, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(ep->Y, y5, ep->Y, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(z5, t2, z5, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(ep->Z, ep->Z, z5, ec->p, ec->n, ec->size);
}
