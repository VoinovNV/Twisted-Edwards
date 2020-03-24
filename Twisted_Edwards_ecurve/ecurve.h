#pragma once
#include <stdio.h>
#include <ak_mpzn.h>
#include <ak_curves.h>

typedef struct ecurve* ak_ecurve;
typedef struct epoint* ak_epoint;

struct epoint {
	ak_uint64 X[ak_mpzn512_size];
	ak_uint64 Y[ak_mpzn512_size];
	ak_uint64 Z[ak_mpzn512_size];
};
struct ecurve {
	ak_uint64 e[ak_mpzn512_size];
	ak_uint64 d[ak_mpzn512_size];
	ak_uint64 p[ak_mpzn512_size];
	ak_uint64 q[ak_mpzn512_size];
	struct epoint point;
	ak_uint32 size;
	ak_uint64 n;
	ak_uint64 r2q[ak_mpzn512_size];
	ak_uint64 r2[ak_mpzn512_size];
	ak_uint64 nq;
	ak_uint64 S[ak_mpzn512_size];
	ak_uint64 T[ak_mpzn512_size];
	const char *pchar;
};

/* ?????????? ????? ???????? ???????? ???????? */
void ak_epoint_set_as_unit(ak_epoint , ak_ecurve);

/* ?????????? ????? ???????? ?????? ????? */
void ak_epoint_set_epoint(ak_epoint, ak_epoint, ak_ecurve);

/* ?????????? ????? ???????? ?????????? ????? ?????? */
void ak_epoint_set(ak_epoint, ak_ecurve);

/* ?????????? ????? ? ???????? ???? */
void ak_epoint_reduce(ak_epoint , ak_ecurve );

/* ???????? ????? ?? ?????????????? ??????*/
bool_t ak_epoint_is_ok(ak_epoint , ak_ecurve );

/* ???????? ????? ?? ?????? */
void ak_epoint_double(ak_epoint , ak_ecurve );

/* ???????? ???? ????? ?? ?????? */
void ak_epoint_add(ak_epoint , ak_epoint , ak_ecurve );

/* ?????????? ??????? ????? */
void ak_epoint_pow(ak_epoint , ak_epoint , ak_uint64 * , size_t , ak_ecurve );


void ak_curve_change_form(ak_ecurve);
/* ??????? ????? ?? ?????? ???????? ? ????? ?? ?????? ????????????*/
void ak_epoint_to_wpoint(ak_epoint, ak_ecurve, ak_wpoint, ak_wcurve);

void ak_wpoint_to_epoint(ak_wpoint, ak_wcurve, ak_epoint, ak_ecurve);

void ak_epoint_triple(ak_epoint , ak_ecurve );

void ak_epoint_quintuple(ak_epoint , ak_ecurve );

/* ????????? ?????????? ?????? ???????? */
const static struct ecurve gost_3410_2012_512_paramSetC = {
	{ 0x0000000000000239, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000 }, /* e ? ????? ?????????? */
	{ 0x6515a5166d05caf7, 0xae6dc7d439a723d5, 0xdc1c74edcea76671, 0x853a44eed58ae3e5, 0xc84c79f64266472e, 0xa1a4bfeccd0cf540, 0xab899e4c73783aa1, 0xde66ec2f500fc692 }, /* d ? ????? ?????????? */
	{ 0xfffffffffffffdc7, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff }, /* p */
	{ 0x94623cef47f023ed, 0xc8eda9e7a769a126, 0x4c33a9ff5147502c, 0xc98cdba46506ab00, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x3fffffffffffffff }, /* q */
	{
		{ 0x0000000000000012, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000 }, /* px */
		{ 0x600303ee73001a3d, 0x905622c04b2baae7, 0xbf068c5d139732f0, 0x22dd4b650cf789ee, 0x9a56117f7b386695, 0x0fdfb0d01794368d, 0x6b99592b77a01e2a, 0x469af79d1fb1f5e1 }, /* py */
		{ 0x0000000000000001, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000 }  /* pz */
	},
	ak_mpzn512_size, /* size */
	0x58a1f7e6ce0f4c09LL, /* n */
	{ 0xe58fa18ee6ca4eb6, 0xe79280282d956fca, 0xd016086ec2d4f903, 0x542f8f3fa490666a, 0x04f77045db49adc9, 0x314e0a57f445b20e, 0x8910352f3bea2192, 0x394c72054d8503be }, /* r2q */
	{ 0x000000000004f0b1, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000 }, /* r2 */
	0x0ed9d8e0b6624e1bLL, /* nq */
	{ 0xa6ba96ba64be8cb4, 0x94648e0af196370a, 0x88f8e2c48c562663, 0x5eb16ec44a9d4706, 0xcdece1826f666e34, 0x9796d004ccbcc2af, 0x551d986ce321f157, 0x486644f42bfc0e5b }, /* S */
	{ 0xe62e462e6780f788, 0x9d124bf8b44685f8, 0xfa04be27a2713bbd, 0x163460d278ec7b50, 0x76b769a90b110bdd, 0xf0461ffcccd77e35, 0x71ec450cbde95f1a, 0x2511275d3802a118 }, /* T */
	"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDC7"
};
