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
#include <emmintrin.h>
#include <mmintrin.h>
#include <smmintrin.h>
#include "utils.h"
#include "bwt.h"
#include "bwt_aux.h"
#include "kvec.h"

static const uint32_t occ_mask[16] = {
	0xc0000000u, 0xf0000000u, 0xfc000000u, 0xff000000u,
	0xffc00000u, 0xfff00000u, 0xfffc0000u, 0xffff0000u, 
	0xffffc000u, 0xfffff000u, 0xfffffc00u, 0xffffff00u,
	0xffffffc0u, 0xfffffff0u, 0xfffffffcu, 0xffffffffu
};

/*static const uint64_t occ_mask2[32] = {
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
};*/

static const uint64_t n_mask[5] = { 0xfffffffffffffffful, 0xaaaaaaaaaaaaaaaaul, 
		0x5555555555555555ul, 0x0ul, 0xfffffffffffffffful };

static __m128i n_mask_128[3];
static __m64 n_mask_64[9];

void bwt_gen_cnt_table(bwt_t *bwt)
{
	int i, j;
	for (i = 0; i != 256; ++i) {
		uint32_t x = 0;
		for (j = 0; j != 4; ++j)
			x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);
		bwt->cnt_table[i] = x;
	}
	n_mask_64[0] = _mm_cvtsi64x_si64(0xfffffffffffffffful);
	n_mask_64[1] = _mm_cvtsi64x_si64(0xaaaaaaaaaaaaaaaaul);
	n_mask_64[2] = _mm_cvtsi64x_si64(0x5555555555555555ul);
	n_mask_64[3] = _mm_setzero_si64();
	n_mask_64[4] = _mm_cvtsi64x_si64(0xfffffffffffffffful);
	n_mask_64[5] = _mm_cvtsi64x_si64(0x3333333333333333ul);
	n_mask_64[6] = _mm_cvtsi64x_si64(0x0f0f0f0f0f0f0f0ful);
	n_mask_64[7] = _mm_cvtsi64x_si64(0x1555555555555555ul);
	n_mask_64[8] = _mm_cvtsi64x_si64(0x1111111111111111ul);
	n_mask_128[0] = _mm_set1_epi64(n_mask_64[2]);
	n_mask_128[1] = _mm_set1_epi64(n_mask_64[5]);
	n_mask_128[2] = _mm_set1_epi64(n_mask_64[6]);
}

static inline int __occ_aux(uint64_t y)
{
	// reduce nucleotide counting to bits counting
	y = (y >> 1) & y;
	// count the number of 1s in y
	y = (y & 0x1111111111111111ul) + (y >> 2 & 0x1111111111111111ul);
	return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0ful) * 0x101010101010101ul >> 56;
}

#define occ_mask2(n) (0x5555555555555555ul - (0x1555555555555555ul >> ((n&31)<<1)))

//reduce nucleotide to bit counting - after w xor n_mask[c].
#define nucleo_5mask(w) 	(w & (w >> 1) & 0x5555555555555555ul)

//transforms in five stages from bits to number of bits.
//these do not allow addition of several of the previous bit-count stage
#define nucleo_3mask(w)		((w + (w >> 2)) & 0x3333333333333333ul)
#define nucleo_f0mask(w)	((w + (w >> 4)) & 0x0f0f0f0f0f0f0f0ful)

// The combined three of these do the same as w * 0x101010101010101ull >> 56
#define nucleo_ffmask(w)	((w + (w >> 8)) & 0x00ff00ff00ff00fful)
#define nucleo_4fmask(w)	((w + (w >> 16)) & 0x0000ffff0000fffful)
#define nucleo_8fmask(w)	((w + (w >> 32)) & 0x00000000fffffffful)

//for convenience, these macros do several stages in one call.
#define nucleo_upto5mask(p, x, w) ({		\
	w = (p) ^ (x);				\
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

// These combine macros are alternative stages from bits to number of bits,
// these once allow the addition of several of the previous bit-count stage
#define nucleo_combine_3mask(v, w) ({		\
	v = w & 0x3333333333333333ul;		\
	v + ((w ^ v) >> 2);			\
})

#define nucleo_combine_f0mask(v, w) ({		\
	v = w & 0xf0f0f0f0f0f0f0ful;		\
	v + ((w ^ v) >> 4);			\
})

#define nucleo_combine_ffmask(v, w) ({		\
	v = w & 0x00ff00ff00ff00fful;		\
	v + ((w ^ v) >> 8);			\
})

#define nucleo_combine_4fmask(v, w) ({		\
	v = w & 0x0000ffff0000fffful;		\
	v + ((w ^ v) >> 16);			\
})

#define nucleo_combine_8fmask(v, w) ({		\
	v = w & 0x00000000fffffffful;		\
	v + ((w ^ v) >> 32);			\
})

//nucleo_upto5mask
#define _mm_nc_combmask_xxx(x, xor, m, t1, t2) ({		\
	x = _mm_xor_##t1(x, xor);				\
	x = _mm_and_##t1(x, _mm_srli_##t2(x, 1));		\
	_mm_and_##t1(x, m);					\
})

#define _mm_nc_combmask_si64(q, x, m) _mm_nc_combmask_xxx(q, x, m, si64, si64)
#define _mm_nc_combmask_epi128(t, x, m) 				\
	_mm_nc_combmask_xxx(t, x, m, si128, epi64)

// e.g. nucleo_3mask()
#define _mm_nmask_xxx(q, q2, shft, m, t1, t2) ({\
	q2 = _mm_srli_##t2(q, shft);		\
	q = _mm_add_##t2(q, q2);		\
	_mm_and_##t1(q, m);			\
})

#define _mm_nmask_si64(q, q2, shft, m)		\
	_mm_nmask_xxx(q, q2, shft, m, si64, si64)

#define _mm_nmask_epi128(q, q2, shft, m)	\
	_mm_nmask_xxx(q, q2, shft, m, si128, epi64)

// e.g. als nucleo_combine_3mask()
#define _mm_cn_mask_xxx(x, y, m, shft, t1, t2) ({	\
	x = _mm_and_##t1(y, m);				\
	y = _mm_xor_##t1(y, x);				\
	y = _mm_srli_##t2(y, shft);			\
	_mm_add_##t2(x, y);				\
})

#define _mm_cn_mask_si64(x, q, m, shft) _mm_cn_mask_xxx(x, q, m, shft, si64, si64)
#define _mm_cn_mask_epi128(t2, t, m, shft) _mm_cn_mask_xxx(t2, t, m, shft, si128, epi64)

#define _mm_sum2si128_si64(t)			\
	_mm_add_si64(_mm_movepi64_pi64(t), _mm_cvtsi64x_si64(_mm_extract_epi64(t, 1)))

static inline uint64_t bwt_occ(bwtint_t k, __m64 x, const uint32_t *const p)
{
	__m128i t, t2;
	t2 = t = _mm_set1_epi64(x);
	switch (k&0x60) {
		case 0x60:
			x = _mm_set_pi32(p[-6], p[-5]);
		case 0x40:
			t = _mm_set_epi64(x, _mm_set_pi32(p[-4], p[-3]));
		case 0x20: x = _mm_set_pi32(p[-2], p[-1]);
	}
	t = _mm_nc_combmask_epi128(t, t2, n_mask_128[0]);

	t2 = _mm_xor_si128(_mm_set_epi64(x, _mm_set_pi32(p[0], p[1])), t2);
	t2 = _mm_and_si128(t2, _mm_srli_epi64(t2, 1));

	x = _mm_srli_si64(n_mask_64[7], ((k&31)<<1));
	x = _mm_sub_si64(n_mask_64[2], x);
	t2 = _mm_and_si128(t2, _mm_set_epi64(n_mask_64[2], x));

	t = _mm_add_epi64(t, t2);

	t2 = _mm_srli_epi64(t, 2);
	t2 = _mm_and_si128(t2, n_mask_128[1]);
	t = _mm_and_si128(t, n_mask_128[1]);
	t = _mm_add_epi64(t, t2);

	t2 = _mm_srli_epi64(t, 4);
	t2 = _mm_and_si128(t2, n_mask_128[2]);
	t = _mm_and_si128(t, n_mask_128[2]);
	t = _mm_add_epi64(t, t2);

	return _mm_extract_epi64(t, 1) + _mm_extract_epi64(t, 0);
}

static inline bwtint_t bwt_invPsi(const bwt_t *bwt, bwtint_t isa)
{
	if (likely(isa != bwt->primary)) {
		bwtint_t c, _i;
		_i = (isa < bwt->primary) ? isa : isa - 1;
		c = bwt_B0(bwt, _i);
		if (likely(isa < bwt->seq_len)) {
			const uint32_t *p;
			isa = bwt->L2[c] + ((const bwtint_t *)(p = bwt_occ_intv(bwt, _i)))[c];
			p += sizeof(bwtint_t) + ((_i&0x60)>>4);
			c = bwt_occ(_i, n_mask_64[c], p);
			isa += c * 0x101010101010101ul >> 56;
		} else {
			isa = (isa == bwt->seq_len ? bwt->L2[c+1] : bwt->L2[c]);
		}
	} else {
		isa = 0;
	}

	return isa;
}

bwtint_t bwt_sa(const bwt_t *bwt, bwtint_t k)
{
	bwtint_t mask, sa = 0;
	mask = bwt->sa_intv - 1;
	while(k & mask) {
		++sa;
		k = bwt_invPsi(bwt, k);
	}
	/* without setting bwt->sa[0] = -1, the following line should be
	   changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1) */
	return sa + bwt->sa[k/bwt->sa_intv];
}


// bwt->bwt and bwt->occ must be precalculated
void bwt_cal_sa(bwt_t *bwt, int intv)
{
	bwtint_t isa, sa, i; // S(isa) = sa
	int intv_round = intv;

	kv_roundup32(intv_round);
	xassert(intv_round == intv, "SA sample interval is not a power of 2.");
	xassert(bwt->bwt, "bwt_t::bwt is not initialized.");

	if (bwt->sa) free(bwt->sa);
	bwt->sa_intv = intv;
	bwt->n_sa = (bwt->seq_len + intv) / intv;
	bwt->sa = (bwtint_t*)calloc(bwt->n_sa, sizeof(bwtint_t));
	if (bwt->sa == 0) {
		fprintf(stderr, "[%s] Fail to allocate %.3fMB memory. Abort!\n", __func__, bwt->n_sa * sizeof(bwtint_t) / 1024.0/1024.0);
		abort();
	}
	// calculate SA value
	isa = 0; sa = bwt->seq_len;
	for (i = 0; i < bwt->seq_len; ++i) {
		if (isa % intv == 0)
			bwt->sa[isa/intv] = sa;
		--sa;
		isa = bwt_invPsi(bwt, isa);
	}
	if (isa % intv == 0) bwt->sa[isa/intv] = sa;
	bwt->sa[0] = (bwtint_t)-1; // before this line, bwt->sa[0] = bwt->seq_len
}

// an analogy to bwt_occ() but more efficient, requiring k <= l
inline bwtint_t bwt_2occ_rk(const bwt_t *bwt, bwtint_t k, bwtint_t *l, ubyte_t c)
{
	if (*l >= bwt->primary) {
		if (k > bwt->primary) {
			--k;
		} else if (k == 0) {
			*l = bwt->L2[c+1];
			k = bwt->L2[c] + 1;
			goto out;
		}
		--*l;
	}
	bwtint_t w, y, z;
	const uint32_t *p, *p2;
	bwtint_t n = *l;

	--k;
	*l = ((const bwtint_t *)(p = bwt_occ_intv(bwt, n)))[c] + bwt->L2[c];
	const uint64_t x = (c == 1 ? 0xaaaaaaaaaaaaaaaaul :
		(c == 2 ? 0x5555555555555555ul :
		(c == 3 ? 0x0ul : 0xfffffffffffffffful)));

	p += sizeof(bwtint_t) + ((n&0x60)>>4);
	w = ((uint64_t)p[0]<<32 | p[1]) ^ x;
	y = w & (w >> 1) & occ_mask2(n);
	n = ((n^k)&~31) | ((k&0x60) >> 4);
	w = 0ul;
	z = occ_mask2(k);
	switch (n) {
	case 0x6: w = nucleo_upto5mask(((uint64_t)p[-6]<<32 | p[-5]), x, k);
	case 0x4: w += nucleo_upto5mask(((uint64_t)p[-4]<<32 | p[-3]), x, k);
	case 0x2: w += nucleo_upto5mask(((uint64_t)p[-2]<<32 | p[-1]), x, k);
		w = nucleo_combine_3mask(k, w);
	case 0x0: z &= y;
		if (y == z) {
			k = (bwtint_t)(-1);
			goto out;
		}
		y = nucleo_3mask(y);
		z = nucleo_3mask(z) + w;
		k = *l;
		break;
	case 0x24:w = nucleo_upto5mask(((uint64_t)p[-6]<<32 | p[-5]), x, k);
	case 0x22:w += nucleo_upto5mask(((uint64_t)p[-4]<<32 | p[-3]), x, k);
		w = nucleo_combine_3mask(k, w);
	case 0x20: k = nucleo_upto5mask(((uint64_t)p[-2]<<32 | p[-1]), x, k);
		y += k;
		k &= z;
		z = nucleo_3mask(k);
		y = nucleo_combine_3mask(k, y);
		if (y == z) {
			k = (bwtint_t)(-1);
			goto out;
		}
		z += w;
		k = *l;
		break;
	default:
		k = ((const bwtint_t *)(p2 = bwt_occ_intv(bwt, k)))[c] + bwt->L2[c];
		n &= 0x66;
		p2 += sizeof(bwtint_t) + (n & 0x6); //is really k
		w = ((uint64_t)p2[0]<<32 | p2[1]) ^ x;
		z &= w & (w >> 1);
		w = 0ul;
		switch (n) {
		case 0x06: w = nucleo_upto5mask(((uint64_t)p[-6]<<32 | p[-5]), x, n);
		case 0x26: w += nucleo_upto5mask(((uint64_t)p[-4]<<32 | p[-3]), x, n);
		case 0x46: w += nucleo_upto5mask(((uint64_t)p[-2]<<32 | p[-1]), x, n);
			w = nucleo_combine_3mask(n, w);
		case 0x66: z += nucleo_upto5mask(((uint64_t)p2[-6]<<32 | p2[-5]), x, n);
			z += nucleo_upto5mask(((uint64_t)p2[-4]<<32 | p2[-3]), x, n);
			z = nucleo_combine_3mask(n, z);
			n = nucleo_upto5mask(((uint64_t)p2[-2]<<32 | p2[-1]), x, n);
			z += nucleo_3mask(n);
			break;
		case 0x60: w = nucleo_upto5mask(((uint64_t)p[-6]<<32 | p[-5]), x, n);
		case 0x40: w += nucleo_upto5mask(((uint64_t)p[-4]<<32 | p[-3]), x, n);
		case 0x20: w += nucleo_upto5mask(((uint64_t)p[-2]<<32 | p[-1]), x, n);
			w = nucleo_combine_3mask(n, w);
		case 0x00: z = nucleo_3mask(z);
			break;
		default:
			switch (n) {
			case 0x24: w = nucleo_upto5mask(((uint64_t)p[-6]<<32 | p[-5]), x, n);
			case 0x04: w += nucleo_upto5mask(((uint64_t)p[-4]<<32 | p[-3]), x, n);
			case 0x64: w += nucleo_upto5mask(((uint64_t)p[-2]<<32 | p[-1]), x, n);
				w = nucleo_combine_3mask(n, w);
			case 0x44: z += nucleo_upto5mask(((uint64_t)p2[-4]<<32 | p2[-3]), x, n);
				break;
			case 0x42: w = nucleo_upto5mask(((uint64_t)p[-6]<<32 | p[-5]), x, n);
			case 0x62: w += nucleo_upto5mask(((uint64_t)p[-4]<<32 | p[-3]), x, n);
			case 0x02: w += nucleo_upto5mask(((uint64_t)p[-2]<<32 | p[-1]), x, n);
				w = nucleo_combine_3mask(n, w);
			}
			z += nucleo_upto5mask(((uint64_t)p2[-2]<<32 | p2[-1]), x, n);
			z = nucleo_combine_3mask(n, z);
		}
		y = nucleo_3mask(y);
	}
	y += w;
	z = nucleo_combine_f0mask(w, z);
	y = nucleo_combine_f0mask(w, y);
	/*y = nucleo_ffmask(y);
	z = nucleo_ffmask(z);
	z = nucleo_4fmask(z);
	y = nucleo_4fmask(y);
	*l += nucleo_8fmask(y);
	k += nucleo_8fmask(z) + 1;*/
	k += (z * 0x101010101010101ul >> 56) + 1;
	*l += y * 0x101010101010101ul >> 56;
out:
	return k;
}

inline bwtint_t bwt_2occ(const bwt_t *bwt, bwtint_t k, bwtint_t *l, ubyte_t c)
{
  bwtint_t ok1, ol1, ok2, ol2, ok3, ol3;
  bwtint_t cntk[4], cntl[4];

  ok1 = ok2 = ok3 = 0;
  ol1 = ol2 = ol3 = 0;
  
  bwt_aux_2occ(bwt, k, *l, c, &ok1, &ol1);
  // NB: we must add the counts, as this is added in the bwt_2occ_rk function
  ok1 += bwt->L2[c] + 1;
  ol1 += bwt->L2[c];

  ol2 = *l;
  ok2 = bwt_2occ_rk(bwt, k, &ol2, c);

  bwt_2occ4(bwt, k, *l, cntk, cntl);
  // NB: we must add the counts, as this is added in the bwt_2occ_rk function
  ok3 = cntk[c] + bwt->L2[c] + 1;
  ol3 = cntl[c] + bwt->L2[c];

  if(ok1 != ok2 || ol1 != ol2 || ok1 != ok3 || ol1 != ol3) {
      fprintf(stderr, "\nERROR:\nc=%u k=%lld *l=%lld\n",
              c, k, *l); 
      fprintf(stderr, "ok1=%lld ok2=%lld ok3=%lld\n",
              ok1, ok2, ok3);
      fprintf(stderr, "ol1=%lld ol2=%lld ol3=%lld\n", 
              ol1, ol2, ol3);
      exit(1);
  }

  *l = ol2;
  return ok2;
  
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
	memcpy(cnt, p, 4 * sizeof(bwtint_t));
	p += sizeof(bwtint_t);
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
	_k = (k >= bwt->primary)? k-1 : k;
	_l = (l >= bwt->primary)? l-1 : l;
	if (_l/OCC_INTERVAL != _k/OCC_INTERVAL || k == (bwtint_t)(-1) ||  unlikely(l == (bwtint_t)(-1))) {
		bwt_occ4(bwt, k, cntk);
		bwt_occ4(bwt, l, cntl);
	} else {
		bwtint_t i, j, x, y;
		uint32_t *p;
		if (k >= bwt->primary) --k; // because $ is not in bwt
		if (l >= bwt->primary) --l;
		p = bwt_occ_intv(bwt, k);
		memcpy(cntk, p, 4 * sizeof(bwtint_t));
		p += sizeof(bwtint_t);
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
		memcpy(cntl, cntk, 4 * sizeof(bwtint_t));
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

/*********************
 * Bidirectional BWT *
 *********************/

void bwt_extend(const bwt_t *bwt, const bwtintv_t *ik, bwtintv_t ok[4], int is_back)
{
	bwtint_t tk[4], tl[4];
	int i;
	bwt_2occ4(bwt, ik->x[!is_back] - 1, ik->x[!is_back] - 1 + ik->x[2], tk, tl);
	for (i = 0; i != 4; ++i) {
		ok[i].x[!is_back] = bwt->L2[i] + 1 + tk[i];
		ok[i].x[2] = tl[i] - tk[i];
	}
	ok[3].x[is_back] = ik->x[is_back] + (ik->x[!is_back] <= bwt->primary && ik->x[!is_back] + ik->x[2] - 1 >= bwt->primary);
	ok[2].x[is_back] = ok[3].x[is_back] + ok[3].x[2];
	ok[1].x[is_back] = ok[2].x[is_back] + ok[2].x[2];
	ok[0].x[is_back] = ok[1].x[is_back] + ok[1].x[2];
}

static void bwt_reverse_intvs(bwtintv_v *p)
{
	if (p->n > 1) {
		int j;
		for (j = 0; j < p->n>>1; ++j) {
			bwtintv_t tmp = p->a[p->n - 1 - j];
			p->a[p->n - 1 - j] = p->a[j];
			p->a[j] = tmp;
		}
	}
}

int bwt_smem1(const bwt_t *bwt, int len, const uint8_t *q, int x, bwtintv_v *mem, bwtintv_v *tmpvec[2])
{
	int i, j, c, ret;
	bwtintv_t ik, ok[4];
	bwtintv_v a[2], *prev, *curr, *swap;

	mem->n = 0;
	if (q[x] > 3) return x + 1;
	kv_init(a[0]); kv_init(a[1]);
	prev = tmpvec[0]? tmpvec[0] : &a[0];
	curr = tmpvec[1]? tmpvec[1] : &a[1];
	bwt_set_intv(bwt, q[x], ik);
	ik.info = x + 1;

	for (i = x + 1, curr->n = 0; i < len; ++i) { // forward search
		if (q[i] < 4) {
			c = 3 - q[i];
			bwt_extend(bwt, &ik, ok, 0);
			if (ok[c].x[2] != ik.x[2]) // change of the interval size
				kv_push(bwtintv_t, *curr, ik);
			if (ok[c].x[2] == 0) break; // cannot be extended
			ik = ok[c]; ik.info = i + 1;
		} else { // an ambiguous base
			kv_push(bwtintv_t, *curr, ik);
			break; // cannot be extended; in this case, i<len always stands
		}
	}
	if (i == len) kv_push(bwtintv_t, *curr, ik); // push the last interval if we reach the end
	bwt_reverse_intvs(curr); // s.t. smaller intervals visited first
	ret = curr->a[0].info; // this will be the returned value
	swap = curr; curr = prev; prev = swap;

	for (i = x - 1; i >= -1; --i) { // backward search for MEMs
		if (q[i] > 3) break;
		c = i < 0? 0 : q[i];
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			bwtintv_t *p = &prev->a[j];
			bwt_extend(bwt, p, ok, 1);
			if (ok[c].x[2] == 0 || i == -1) { // keep the hit if reaching the beginning or not extended further
				if (curr->n == 0) { // curr->n to make sure there is no longer matches
					if (mem->n == 0 || i + 1 < mem->a[mem->n-1].info>>32) { // skip contained matches
						ik = *p; ik.info |= (uint64_t)(i + 1)<<32;
						kv_push(bwtintv_t, *mem, ik);
					}
				} // otherwise the match is contained in another longer match
			}
			if (ok[c].x[2] && (curr->n == 0 || ok[c].x[2] != curr->a[curr->n-1].x[2])) {
				ok[c].info = p->info;
				kv_push(bwtintv_t, *curr, ok[c]);
			}
		}
		if (curr->n == 0) break;
		swap = curr; curr = prev; prev = swap;
	}
	bwt_reverse_intvs(mem); // s.t. sorted by the start coordinate

	if (tmpvec[0] == 0) free(a[0].a);
	if (tmpvec[1] == 0) free(a[1].a);
	return ret;
}
