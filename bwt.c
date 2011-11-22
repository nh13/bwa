/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include "utils.h"
#include "bwt.h"

static const uint32_t occ_mask[16] = {
	0xc0000000u, 0xf0000000u, 0xfc000000u, 0xff000000u,
	0xffc00000u, 0xfff00000u, 0xfffc0000u, 0xffff0000u, 
	0xffffc000u, 0xfffff000u, 0xfffffc00u, 0xffffff00u,
	0xffffffc0u, 0xfffffff0u, 0xfffffffcu, 0xffffffffu
};

static const uint64_t occ_mask2[32] = {
	0x40000000ul, 	0x50000000ul, 	0x54000000ul,
	0x55000000ul, 	0x55400000ul, 	0x55500000ul,
	0x55540000ul, 	0x55550000ul, 	0x55554000ul,
	0x55555000ul, 	0x55555400ul, 	0x55555500ul,
	0x55555540ul, 	0x55555550ul, 	0x55555554ul,
	0x55555555ul, 	0x4000000055555555ul, 	0x5000000055555555ul,
	0x5400000055555555ul, 	0x5500000055555555ul, 	0x5540000055555555ul,
	0x5550000055555555ul, 	0x5554000055555555ul, 	0x5555000055555555ul,
	0x5555400055555555ul, 	0x5555500055555555ul, 	0x5555540055555555ul,
	0x5555550055555555ul, 	0x5555554055555555ul, 	0x5555555055555555ul,
	0x5555555455555555ul, 0x5555555555555555ul
};

static const uint64_t n_mask[5] = { 0xfffffffffffffffful, 0xaaaaaaaaaaaaaaaaul, 
		0x5555555555555555ul, 0x0ul, 0xfffffffffffffffful };

void bwt_gen_cnt_table(bwt_t *bwt)
{
	int i, j;
	for (i = 0; i != 256; ++i) {
		uint32_t x = 0;
		for (j = 0; j != 4; ++j)
			x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);
		bwt->cnt_table[i] = x;
	}
}

static inline int __occ_aux(uint64_t y)
{
	// reduce nucleotide counting to bits counting
	y = (y >> 1) & y;
	// count the number of 1s in y
	y = (y & 0x1111111111111111ul) + (y >> 2 & 0x1111111111111111ul);
	return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0ful) * 0x101010101010101ul >> 56;
}

#define nucleo_5mask(w) 	(w & (w >> 1) & 0x5555555555555555ul)
#define nucleo_3mask(w)		((w + (w >> 2)) & 0x3333333333333333ul)
#define nucleo_f0mask(w)	((w + (w >> 4)) & 0x0f0f0f0f0f0f0f0ful)

#define nucleo_upto5mask(p, x, w) ({		\
	w = *(p) ^ (x);				\
	nucleo_5mask(w);			\
})

#define nucleo_upto3mask(p, x, w) ({		\
	w = nucleo_upto5mask(p, x, w);		\
	nucleo_3mask(w);			\
})

#define nucleo_uptof0mask(p, x, w) ({		\
	w = nucleo_upto3mask(p, x, w);		\
	nucleo_f0mask(w);			\
})


#define nucleo_pre_combine_3mask(v, w) ({	\
	v = w;					\
	w &= 0x3333333333333333ul;		\
})

// doesn't work?
#define nucleo_pre_combine_f0mask(v, w) ({	\
	v = w;					\
	w &= 0x0f0f0f0f0f0f0f0ful;		\
	w ^ v;					\
})

#define nucleo_combine_3mask(v, w) ({		\
	v = w & 0x3333333333333333ul;		\
	v + ((w ^ v) >> 2);			\
})

#define nucleo_combine_f0mask(v, w) ({		\
	v = w & 0xf0f0f0f0f0f0f0ful;		\
	v + ((w ^ v) >> 4);			\
})

static inline uint64_t bwt_occ(const bwtint_t k, uint64_t x, const uint64_t *const p)
{
	uint64_t w = 0ul;
	uint64_t y = *p ^ x;
	y &= (y >> 1) & occ_mask2[k&31];
	switch (k&0x60) {
		case 0x60: y += nucleo_upto5mask(p - 3, x, w);
			y += nucleo_upto5mask(p - 2, x, w);
			y = nucleo_combine_3mask(w, y);
			y += nucleo_upto3mask(p - 1, x, w);
			return nucleo_combine_f0mask(w, y);
		case 0x00: y = nucleo_3mask(y);
			break;
		case 0x40: y += nucleo_upto5mask(p - 2, x, w);
		case 0x20: y += nucleo_upto5mask(p - 1, x, w);
			y = nucleo_combine_3mask(w, y);
	}
	return nucleo_f0mask(y);
}

static inline bwtint_t cal_isa(const bwt_t *bwt, bwtint_t isa)
{
	if (likely(isa != bwt->primary)) {
		bwtint_t c, _i, so;
		_i = (isa < bwt->primary) ? isa : isa - 1;
		so = _i/OCC_INTERVAL;
		c = bwt_B0(bwt, _i, so);
		if (likely(isa < bwt->seq_len)) {
			uint64_t w;
			const uint64_t *p = (const uint64_t *)bwt->bwt + so * 6;
			w = n_mask[c];
			isa = bwt->L2[c] + ((uint32_t *)p)[c];
			p += 2 + ((_i&0x60)>>5);
			w = bwt_occ(_i, w, p);
			isa += w * 0x101010101010101ul >> 56;
		} else {
			isa = (isa == bwt->seq_len ? bwt->L2[c+1] : bwt->L2[c]);
		}
	} else {
		isa = 0;
	}

	return isa;
}

static inline bwtint_t cal_isa_PleSl(const bwt_t *bwt, bwtint_t isa)
{
	uint64_t w;
	const uint64_t *p;
	bwtint_t c;
	c = isa/OCC_INTERVAL; //unaffected by isa decr.
	p = (const uint64_t *)bwt->bwt + c * 6;
	if (likely(isa != bwt->primary)) {
		if (isa > bwt->primary) {
			if (unlikely(isa > bwt->seq_len))
				return 0;
			--isa;

		}
		c = bwt_B0(bwt, isa, c);
		w = n_mask[c];
		c = bwt->L2[c] + ((uint32_t *)p)[c];
		p += 2 + ((isa&0x60)>>5);
		w = bwt_occ(isa, w, p);
		isa = c + (w * 0x101010101010101ul >> 56);
		 //only 0x1f1f... part in _i&31?
		 //can we reduce this since we only really have to count the first bits?
	} else {
		c = bwt_B0(bwt, isa, c);
		if (isa == bwt->seq_len)
			++c;
		isa = bwt->L2[c];
	}
	return isa;
}

bwtint_t bwt_sa(const bwt_t *bwt, bwtint_t k)
{
	bwtint_t m, sa = 0;
	m = bwt->sa_intv - 1;
	if (likely(!((m+1) & m) && bwt->primary <= bwt->seq_len)) {
		// not power of 2 before decrement
		while (k & m) {
			++sa;
			k = cal_isa_PleSl(bwt, k);
		}
	} else {
		bwtint_t add = m;
		m |= m>>1;
		m |= m>>2;
		m |= m>>4;
		m |= m>>8;
		m |= m>>16;
		add ^= m;
		while((k + add) & m) {
			++sa;
			k = cal_isa(bwt, k);
		}
	}
	//k += (z * 0x101010101010101ul >> 56);
	/* without setting bwt->sa[0] = -1, the following line should be
	   changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1) */
	return sa + bwt->sa[k/bwt->sa_intv];
}


// bwt->bwt and bwt->occ must be precalculated
void bwt_cal_sa(bwt_t *bwt, int intv)
{
	bwtint_t isa, sa, i; // S(isa) = sa

	xassert(bwt->bwt, "bwt_t::bwt is not initialized.");

	if (bwt->sa) free(bwt->sa);
	bwt->sa_intv = intv;
	bwt->n_sa = (bwt->seq_len + intv) / intv;
	bwt->sa = (bwtint_t*)calloc(bwt->n_sa, sizeof(bwtint_t));
	// calculate SA value
	isa = 0; sa = bwt->seq_len;
	for (i = 0; i < bwt->seq_len; ++i) {
		if (isa % intv == 0)
			bwt->sa[isa/intv] = sa;
		--sa;
		isa = cal_isa(bwt, isa);
	}
	if (isa % intv == 0) bwt->sa[isa/intv] = sa;
	bwt->sa[0] = (bwtint_t)-1; // before this line, bwt->sa[0] = bwt->seq_len
}

// an analogy to bwt_occ() but more efficient, requiring k <= l
inline bwtint_t bwt_2occ(const bwt_t *bwt, bwtint_t k, bwtint_t *l, ubyte_t c)
{
	if (*l >= bwt->primary) {
		if (k > bwt->primary) {
			--k;
		} else if (k == 0) { // z of previous == 1ul
			*l = bwt->L2[c+1];
			k = bwt->L2[c] + 1;
			goto out;
		}
		--*l;
	}
	register uint64_t v;
	uint64_t w, y, z;
	const uint64_t *p, *p2;
	bwtint_t n = *l;

	--k;
	p = (const uint64_t *)bwt->bwt + n/OCC_INTERVAL * 6;
	*l = ((uint32_t *)p)[c] + bwt->L2[c];
	uint64_t x = (c == 1 ? 0xaaaaaaaaaaaaaaaaul :
		(c == 2 ? 0x5555555555555555ul :
		(c == 3 ? 0x0ul : 0xfffffffffffffffful)));

	w = 0ul;
	p += 2 + ((n&0x60)>>5);
	v = *p ^ x;
	y = v & (v >> 1) & occ_mask2[n&31];
	z = occ_mask2[k&31];
	switch (((n&~31) - (k&~31)) | ((n&0x60) >> 5)) {
	case 0x3: w = nucleo_upto5mask(p - 3, x, v);
	case 0x2: w += nucleo_upto5mask(p - 2, x, v);
	case 0x1: w += nucleo_upto5mask(p - 1, x, v);
		v ^= nucleo_pre_combine_3mask(v, w);
		w += (v >> 2);
	case 0x0: z &= y;
		if (y == z) {
			k = (bwtint_t)(-1);
			goto out;
		}
		y = nucleo_3mask(y);
		z = nucleo_3mask(z) + w;
		k = *l;
		break;
	case 0x23:w = nucleo_upto5mask(p - 3, x, v);
	case 0x22:w += nucleo_upto5mask(p - 2, x, v);
		v ^= nucleo_pre_combine_3mask(v, w);
		w += (v >> 2);
	case 0x21: v = nucleo_upto5mask(p - 1, x, v);
		y += v;
		v &= z;
		z = nucleo_3mask(v);
		v ^= nucleo_pre_combine_3mask(v, y);
		y += (v >> 2);
		if (y == z) {
			k = (bwtint_t)(-1);
			goto out;
		}
		z += w;
		k = *l;
		break;
	default:
		p2 = (const uint64_t *)bwt->bwt + k/OCC_INTERVAL * 6;
		n = (n & 0x60) | ((k & 0x60) >> 5);
		k = ((uint32_t *)p2)[c] + bwt->L2[c];
		p2 += 2 + (n & 0x3); //is really k

		v = *p2 ^ x;
		z &= v & (v >> 1);
		switch (n) {
		case 0x63: w = nucleo_upto5mask(p - 3, x, v);
		case 0x43: w += nucleo_upto5mask(p - 2, x, v);
		case 0x23: w += nucleo_upto5mask(p - 1, x, v);
			v ^= nucleo_pre_combine_3mask(v, w);
			w += (v >> 2);
		case 0x03: z += nucleo_upto5mask(p2 - 3, x, v);
			z += nucleo_upto5mask(p2 - 2, x, v);
			x = nucleo_upto5mask(p2 - 1, x, x);
			v ^= nucleo_pre_combine_3mask(v, z);
			z += (v >> 2) + nucleo_3mask(x);
			break;
		case 0x60: w = nucleo_upto5mask(p - 3, x, v);
		case 0x40: w += nucleo_upto5mask(p - 2, x, v);
		case 0x20: w += nucleo_upto5mask(p - 1, x, v);
			v ^= nucleo_pre_combine_3mask(v, w);
			w += (v >> 2);
		case 0x00: z = nucleo_3mask(z);
			break;
		default:
			switch (n) {
			case 0x62: w = nucleo_upto5mask(p - 3, x, v);
			case 0x42: w += nucleo_upto5mask(p - 2, x, v);
			case 0x22: w += nucleo_upto5mask(p - 1, x, v);
				v ^= nucleo_pre_combine_3mask(v, w);
				w += (v >> 2);
			case 0x02: z += nucleo_upto5mask(p2 - 2, x, v);
				break;
			case 0x61: w = nucleo_upto5mask(p - 3, x, v);
			case 0x41: w += nucleo_upto5mask(p - 2, x, v);
			case 0x21: w += nucleo_upto5mask(p - 1, x, v);
				v ^= nucleo_pre_combine_3mask(v, w);
				w += (v >> 2);
			}
			z += nucleo_upto5mask(p2 - 1, x, v);
			v ^= nucleo_pre_combine_3mask(v, z);
			z += (v >> 2);
		}
		y = nucleo_3mask(y);
	}
	y += w;
	//w = nucleo_pre_combine_f0mask(w, z);
	//z += (w >> 2); //werkt niet ?!
	z = nucleo_combine_f0mask(w, z);
	y = nucleo_combine_f0mask(v, y);
	k += (z * 0x101010101010101ul >> 56) + 1;
	*l += y * 0x101010101010101ul >> 56;
out:
	return k;
}

#define __occ_aux4(bwt, b)											\
	((bwt)->cnt_table[(b)&0xff] + (bwt)->cnt_table[(b)>>8&0xff]		\
	 + (bwt)->cnt_table[(b)>>16&0xff] + (bwt)->cnt_table[(b)>>24])

inline void bwt_occ4(const bwt_t *bwt, bwtint_t k, bwtint_t cnt[4])
{
	bwtint_t l, j, x;
	uint32_t *p;
	if (k == (bwtint_t)(-1)) {
		memset(cnt, 0, 4 * sizeof(bwtint_t));
		return;
	}
	if (k >= bwt->primary) --k; // because $ is not in bwt
	p = bwt_occ_intv(bwt, k);
	memcpy(cnt, p, 16);
	p += 4;
	j = k >> 4 << 4;
	for (l = k / OCC_INTERVAL * OCC_INTERVAL, x = 0; l < j; l += 16, ++p)
		x += __occ_aux4(bwt, *p);
	x += __occ_aux4(bwt, *p & occ_mask[k&15]) - (~k&15);
	cnt[0] += x&0xff; cnt[1] += x>>8&0xff; cnt[2] += x>>16&0xff; cnt[3] += x>>24;
}

// an analogy to bwt_occ4() but more efficient, requiring k <= l
inline void bwt_2occ4(const bwt_t *bwt, bwtint_t k, bwtint_t l, bwtint_t cntk[4], bwtint_t cntl[4])
{
	bwtint_t _k, _l;
	if (k == l) {
		bwt_occ4(bwt, k, cntk);
		memcpy(cntl, cntk, 4 * sizeof(bwtint_t));
		return;
	}
	_k = (k >= bwt->primary)? k-1 : k;
	_l = (l >= bwt->primary)? l-1 : l;
	if (_l/OCC_INTERVAL != _k/OCC_INTERVAL || k == (bwtint_t)(-1) ||  unlikely(l == (bwtint_t)(-1))) {
		bwt_occ4(bwt, k, cntk);
		bwt_occ4(bwt, l, cntl);
	} else {
		bwtint_t i, j, x, y;
		uint32_t *p;
		int cl[4];
		if (k >= bwt->primary) --k; // because $ is not in bwt
		if (l >= bwt->primary) --l;
		cl[0] = cl[1] = cl[2] = cl[3] = 0;
		p = bwt_occ_intv(bwt, k);
		memcpy(cntk, p, 4 * sizeof(bwtint_t));
		p += 4;
		// prepare cntk[]
		j = k >> 4 << 4;
		for (i = k / OCC_INTERVAL * OCC_INTERVAL, x = 0; i < j; i += 16, ++p)
			x += __occ_aux4(bwt, *p);
		y = x;
		x += __occ_aux4(bwt, *p & occ_mask[k&15]) - (~k&15);
		// calculate cntl[] and finalize cntk[]
		j = l >> 4 << 4;
		for (; i < j; i += 16, ++p)
			y += __occ_aux4(bwt, *p);
		y += __occ_aux4(bwt, *p & occ_mask[l&15]) - (~l&15);
		memcpy(cntl, cntk, 16);
		cntk[0] += x&0xff; cntk[1] += x>>8&0xff; cntk[2] += x>>16&0xff; cntk[3] += x>>24;
		cntl[0] += y&0xff; cntl[1] += y>>8&0xff; cntl[2] += y>>16&0xff; cntl[3] += y>>24;
	}
}

int bwt_match_exact(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *sa_begin, bwtint_t *sa_end)
{
	bwtint_t k, l;
	int i;
	k = 0; l = bwt->seq_len;
	for (i = len - 1; i >= 0; --i) {
		ubyte_t c = str[i];
		if (c > 3 || (k = bwt_2occ(bwt, k, &l, c)) == (bwtint_t)(-1))
			return 0; // no match
	}
	if (sa_begin) *sa_begin = k;
	if (sa_end)   *sa_end = l;
	return l - k + 1;
}

int bwt_match_exact_alt(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *k0, bwtint_t *l0)
{
	int i;
	bwtint_t k, l;
	k = *k0; l = *l0;
	for (i = len - 1; i >= 0; --i) {
		ubyte_t c = str[i];
		if (unlikely(c > 3) || (k = bwt_2occ(bwt, k, &l, c)) == (bwtint_t)(-1))
			return 0; // no match
	}
	*k0 = k; *l0 = l;
	return l - k + 1;
}
