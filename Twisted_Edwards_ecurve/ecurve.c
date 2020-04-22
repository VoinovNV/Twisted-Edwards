#include "ecurve.h"
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
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
    ak_mpznmax x5,y5,z5,t1,t2,t_;
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
/////////////////////////////////////////////////////////////////////////////////////////////////////
void ak_epoint_pow_binary( ak_epoint wq, ak_epoint wp, ak_uint64 *k, size_t size,
                          ak_ecurve ec )
{
    ak_uint64 uk = 0;
    long long int i, j;
    struct epoint Q, R;
    ak_epoint_set_as_unit(&Q, ec);
    ak_epoint_set_epoint( &R, wp, ec );
    for( i = size-1; i >= 0; i-- ) {
        uk = k[i];
        for (j = 0; j < 64; j++) {
            ak_epoint_double( &Q, ec );
            if (uk & 0x8000000000000000LL) ak_epoint_add( &Q, &R, ec );
            uk <<= 1;
        }
    }
    ak_epoint_set_epoint( wq, &Q, ec );
}
ak_uint32 ak_mpzn_rem_uint32( ak_uint64 *x, const size_t size, ak_uint32 p )
{
    size_t i;
    ak_uint64 t, r1, r = 1, sum = x[0]%p;

    if( !p ) return ak_error_message( ak_error_invalid_value, __func__, "divide by zero" );

    /* вычисляем r[i] = 2^{64i} mod (p) */
    r1 = 9223372036854775808ull % p; r1 = ( 2*r1 )%p;

    /* вычисляем остаток от деления */
    for( i = 1; i < size; i++ ) {
        r *= r1; r %= p;
        t = x[i]%p; t *= r; t %= p; sum += t;
    }
    /* приведение после суммирования дает корректный результат только для небольших значений size
     в произвольном случае здесь возникнет ошибка переполнения */
    return sum%p;
}

ak_int32 ak_mods(ak_uint64* n,ak_uint32 l,ak_uint32 w,ak_ecurve ec){
    /*ak_int32 a=(pow(l,w));
    if (n%a>=a/l) return n%a-a;
    return n%a;*/
//n mods l^w = (n mod l^w) - l^w , n mod l^w >=l^(w-1)
// else: n mods l^w = n mod l^w
    ak_uint32 a=(pow(l,w-1));
    ak_uint32 p=a*l;
    ak_uint32 res=ak_mpzn_rem_uint32(n, ec->size, p); //result = n mod l^w
    if(res>=a) return res-p;
    return res;
}
void ak_mpzn_div2(ak_uint64* n, size_t size){
    for(unsigned i=0; i<size-1; i++){
        n[i]>>=1;
        if(n[i+1]&1) n[i]|=0x8000000000000000LL;
    }
    n[size-1]>>=1;
}
ak_int32 ak_n_to_NAF2(ak_uint64 *k, ak_int32* res, size_t size){
    /*
     * Atomicity Improvement for Elliptic Curve Scalar Multiplication
     */
    ak_int32 i=0;
    ak_mpznmax E,one=ak_mpznmax_one;
    ak_mpzn_set(E,k,size);
    while (!ak_mpzn_cmp_ui(E,size,0)){
        if(E[0]&1) {
            res[i]=2-(E[0]&3);
            if (res[i]==1) E[0]-=1;
            else (ak_mpzn_add(E,E,one, size));
        }
        else res[i]=0;
        ak_mpzn_div2(E,size);
        i++;
    }
    return i;
}
void ak_invert_ep(ak_epoint resp,ak_epoint ep, ak_ecurve ec){
    ak_mpzn_sub(resp->X,ec->p,ep->X,ec->size);
    memcpy(resp->Y, ep->Y, ec->size * sizeof(ak_uint64));
    memcpy(resp->Z, ep->Z, ec->size * sizeof(ak_uint64));
}
void ak_epoint_pow_NAF( ak_epoint wq, ak_epoint wp, ak_uint64 *k, size_t size,
                       ak_ecurve ec)
{
    ak_int32 kNAF[512];
    ak_int32 i=ak_n_to_NAF2(k,kNAF,ec->size);
    struct epoint Q, R,R_;
    ak_epoint_set_as_unit( &Q, ec );
    ak_epoint_set_epoint( &R, wp, ec );
    ak_epoint_set_epoint( &R_, wp, ec );
    ak_invert_ep(&R_,&R,ec);
    for (i=i-1; i>=0; i--) {
        ak_epoint_double(&Q,ec);
        ak_int32 nj=kNAF[i];
        if (nj>0) {
            ak_epoint_add(&Q,&R,ec);
        }
        else if(nj<0){
            ak_epoint_add(&Q,&R_,ec);
        }
    }
    ak_epoint_set_epoint( wq, &Q, ec );
}

ak_int32 ak_n_to_NAF_powof2(ak_uint64 *k,ak_int32 w, ak_int32* res, size_t size){
    ak_int32 i=0,mods=pow(2,w);
    ak_mpznmax E,one=ak_mpznmax_one;
    ak_mpzn_set(E,k,size);
    while (!ak_mpzn_cmp_ui(E,size,0)){
        if(E[0]&1){
            res[i]=(E[0]&(mods-1));
            if(res[i]>=(mods)>>1) res[i]-=(mods);
            if (res[i]>0){
                ak_mpzn_set_ui(one,size,res[i]);
                ak_mpzn_sub(E,E,one, size);
            }
            else{
                ak_mpzn_set_ui(one,size,-res[i]);
                ak_mpzn_add(E,E,one, size);
            }
        }
        else res[i]=0;
        ak_mpzn_div2(E,size);
        i++;
    }
    return i;
}
void ak_epoint_pow_NAF_powof2( ak_epoint wq, ak_epoint wp, ak_uint64 *k, size_t size,ak_uint32 w,
                       ak_ecurve ec)
{/*
Fixed-Base Comb with Window-Non-Adjacent Form (NAF)
Method for Scalar Multiplication
*/
    ak_mpznmax pows;
    ak_int32 kNAF[512];
    ak_int32 i=ak_n_to_NAF_powof2(k,w,kNAF,ec->size);
    struct epoint Q, R_[32],inv;
    ak_epoint_set_as_unit(&Q, ec );
    ak_epoint_set_epoint(R_, wp, ec );
    ak_epoint_set_epoint(R_+1, wp, ec );
    ak_epoint_set_epoint(R_+2, wp, ec );
    ak_epoint_triple(R_+1, ec);
    if(w>3){
    ak_epoint_quintuple(R_+2, ec);
    for (ak_uint32 j=3;j<w;j++){
        ak_mpzn_set_ui(pows,size,j*2+1);
        ak_epoint_pow_binary(R_+j,R_,pows,size,ec);
    }
    }
    for (i=i-1; i>=0; i--) {
        ak_epoint_double(&Q,ec);
        ak_int32 nj=kNAF[i];
        if (nj>0) {
            ak_epoint_add(&Q,R_+((nj-1)>>1),ec);
        }
        else if(nj<0){
            ak_invert_ep(&inv,R_+((-nj-1)>>1),ec);
            ak_epoint_add(&Q,&inv,ec);
        }
    }
    ak_epoint_set_epoint( wq, &Q, ec );
}


ak_int32 ak_n_to_NAF_L_w(ak_uint64 *k, ak_int32* res, ak_int32 L, ak_int32 w, size_t size,ak_ecurve ec){
    ak_int32 i=0,L_w1=pow(L,w-1);
    ak_int32 L_w=L_w1*L;
    ak_mpznmax E,one,L_inv=ak_mpznmax_zero,u=ak_mpznmax_one,z;
//Находим l^(-1)
    ak_mpzn_set_ui(z, ec->size, 2);
    ak_mpzn_sub(z, ec->p, z, ec->size);
    ak_mpzn_set_ui(L_inv,size,L);
    ak_mpzn_mul_montgomery(L_inv,L_inv,ec->r2,ec->p,ec->n,ec->size);
    ak_mpzn_modpow_montgomery(L_inv, L_inv, z, ec->p, ec->n, ec->size);

    ak_mpzn_set(E,k,size);


    while (ak_mpzn_cmp(E,u,size)>=0){
        if(ak_mpzn_rem_uint32(E,size,L)!=0){
            res[i]=(ak_mpzn_rem_uint32(E,size,L_w));
            if(res[i]>=L_w1) res[i]-=(L_w);
            if (res[i]>0){
                ak_mpzn_set_ui(one,size,res[i]);
                //if(!ak_mpzn_cmp(E,one,size)) break;
                ak_mpzn_sub(E,E,one, size);
            }
            else{
                ak_mpzn_set_ui(one,size,-res[i]);
                ak_mpzn_add(E,E,one, size);
            }
        }
        else res[i]=0;
        ak_mpzn_mul_montgomery(E,E,ec->r2,ec->p,ec->n,ec->size);
        ak_mpzn_mul_montgomery(E,E,L_inv,ec->p,ec->n,ec->size);
        ak_mpzn_mul_montgomery(E,E,u,ec->p,ec->n,ec->size);
        i++;
    }
    return i;
}
void ak_epoint_pow_NAF_powofL( ak_epoint wq, ak_epoint wp, ak_uint64 *k, size_t size, ak_uint32 l,ak_uint32 w,
                              ak_ecurve ec){
    ak_mpznmax pows;
    ak_int32 kNAF[512],LL_=pow(l,w-1);
    ak_int32 len=LL_*(l-1);

    ak_int32 i=ak_n_to_NAF_L_w(k,kNAF,l,w,size,ec);
    struct epoint Q, R_[512],inv;
    ak_epoint_set_as_unit(&Q, ec );
    ak_epoint_set_epoint(R_, wp, ec );
    ak_int32 cr=1;
    for (ak_int32 j=2;j<len;j++){
        if(j%l){ ak_mpzn_set_ui(pows,size,j);
            ak_epoint_pow_binary(R_+cr,R_,pows,size,ec);
            cr++;
        }
    }

    for (i=i-1; i>=0; i--) {
        if(l==2) ak_epoint_double(&Q,ec);
        else{
            if(l==3) ak_epoint_triple(&Q,ec);
            else if(l==5) ak_epoint_quintuple(&Q,ec);
                else{
                    ak_mpzn_set_ui(pows,size,l);
                    ak_epoint_pow_binary(&Q,&Q,pows,size,ec);
                }
        }

        ak_int32 nj=kNAF[i];
        if (nj>0) {
            ak_epoint_add(&Q,R_+(nj-1-(nj/l)),ec);
        }
        else if(nj<0){
            ak_invert_ep(&inv,R_+((-nj)-1-(-nj)/l),ec);
            ak_epoint_add(&Q,&inv,ec);
        }
    }

    ak_epoint_set_epoint( wq, &Q, ec );
}

ak_int32 ak_n_to_Ext_wmb_NAF(ak_uint64 *k, ak_int32* res, ak_int32* bases,
                        ak_uint32* a, ak_uint32* w, ak_int8 len,size_t size,ak_ecurve ec){

    ak_int32 i=0;
    ak_int32 A=1;
    ak_mpznmax E,one,L_inv[10],u=ak_mpznmax_one,z;
    for(int j=0;j<len;j++) A*=pow(a[j],w[j]);

    ak_mpzn_set_ui(z, ec->size, 2);
    ak_mpzn_sub(z, ec->p, z, ec->size);
    for(int j=0;j<len;j++){
        ak_mpzn_set_ui(L_inv[j],size,a[j]);
        ak_mpzn_mul_montgomery(L_inv[j],L_inv[j],ec->r2,ec->p,ec->n,ec->size);
        ak_mpzn_modpow_montgomery(L_inv[j], L_inv[j], z, ec->p, ec->n, ec->size);
    }

    ak_mpzn_set(E,k,size);


    while (ak_mpzn_cmp(E,u,size)>=0){
        int j=0;
        for(;j<len;j++){
            if(ak_mpzn_rem_uint32(E,size,a[j])==0) {res[i]=0; break;}
        }
        if(j==len) {
            res[i]=ak_mpzn_rem_uint32(E,size,A);
            if(res[i]>=A/2) res[i]-=A;
            if (res[i]>0){
                ak_mpzn_set_ui(one,size,res[i]);
                ak_mpzn_sub(E,E,one, size);
            }
            else{
                ak_mpzn_set_ui(one,size,-res[i]);
                ak_mpzn_add(E,E,one, size);
            }
        }

        for(j=0;j<len;j++){
            if(ak_mpzn_rem_uint32(E,size,a[j])==0){
                ak_mpzn_mul_montgomery(E,E,ec->r2,ec->p,ec->n,ec->size);
                ak_mpzn_mul_montgomery(E,E,L_inv[j],ec->p,ec->n,ec->size);
                ak_mpzn_mul_montgomery(E,E,u,ec->p,ec->n,ec->size);
                bases[i]=a[j];
                break;
            }
        }
        i++;
    }

    return i;
}
void ak_epoint_pow_NAF_mbw( ak_epoint wq, ak_epoint wp, ak_uint64 *k, size_t size, ak_uint32* l,ak_uint32* w, ak_uint8 len,
                                   ak_ecurve ec){
    ak_mpznmax pows;
    ak_int32 kNAF[512],bases[512];
    ak_int32 A=1;
    for(int j=0;j<len;j++) A*=pow(l[j],w[j]);
    ak_int32 i=ak_n_to_Ext_wmb_NAF(k,kNAF,bases,l,w,len,size,ec);
    struct epoint Q, R_[13500],inv;
    ak_epoint_set_as_unit(&Q, ec );
    ak_epoint_set_epoint(R_, wp, ec );
    for (ak_int32 j=2;j<A/2;j++){
        ak_mpzn_set_ui(pows,size,j);
        ak_epoint_pow_binary(R_+j-1,R_,pows,size,ec);
    }

    for (i=i-1; i>=0; i--) {

        switch (bases[i]) {
        case 2:
            ak_epoint_double(&Q,ec); break;
        case 3:
            ak_epoint_triple(&Q,ec); break;
        case 5:
            ak_epoint_quintuple(&Q,ec); break;
        default:
            ak_mpzn_set_ui(pows,size,bases[i]);
            ak_epoint_pow_binary(&Q,&Q,pows,size,ec);
            break;
        }

        ak_int32 nj=kNAF[i];
        if (nj>0) {
            ak_epoint_add(&Q,R_+(nj-1),ec);
        }
        else if(nj<0){
            ak_invert_ep(&inv,R_+(-nj-1),ec);
            ak_epoint_add(&Q,&inv,ec);
        }
    }

    ak_epoint_set_epoint( wq, &Q, ec );
}
