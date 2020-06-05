/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * Modifed from gdk_imprints.c
 * Jiading Guo (zingdle@gmail.com)
 */

#include "monetdb_config.h"
#include "gdk.h"
#include "gdk_private.h"
#include "gdk_bindex.h"

#define BINDEX_VERSION		1
#define BINDEX_HEADER_SIZE	3 /* nr of size_t fields in header */

// forward declaration
static void printfv(void *fvi_, BUN bmn);
static void printarea(area_t *area, int tpe);
static void printrv(BUN *rv, BUN cnt);

#define LOCATE_POS(TYPE)	\
	do {	\
		const TYPE *restrict bcol = (TYPE *) Tloc(b, 0);	\
		const TYPE *restrict cmp = (TYPE *) incmp;	\
		area_t *thearea;	\
		BUN first, last, mid;	\
		/* which area */	\
		first = 0;	\
		last = BINDEX_K;	\
		while (last > first) {	\
			mid = first + (last - first) / 2;	\
			if (area[mid].val.TYPE##val < *cmp)	\
				first = mid + 1;	\
			else	\
				last = mid;	\
		}	\
		*al = (first == 0) ? 0 : (first - 1);	\
		first = 0;	\
		last = BINDEX_K;	\
		while (last > first) {	\
			mid = first + (last - first) / 2;	\
			if (area[mid].val.TYPE##val <= *cmp)	\
				first = mid + 1;	\
			else	\
				last = mid;	\
		}	\
		*au = (first == 0) ? 0 : (first - 1);	\
		/* which pos */	\
		thearea = &area[*al];	\
		first = thearea->pos;	\
		last = first + thearea->len;	\
		while (first < last) {	\
			mid = first + (last - first) / 2;	\
			if (bcol[rv[mid]] < *cmp)	\
				first = mid + 1;	\
			else	\
				last = mid;	\
		}	\
		*pl = first;	\
		thearea = &area[*au];	\
		first = thearea->pos;	\
		last = first + thearea->len;	\
		while (first < last) {	\
			mid = first + (last - first) / 2;	\
			if (bcol[rv[mid]] <= *cmp)	\
				first = mid + 1;	\
			else	\
				last = mid;	\
		}	\
		*pu = first;	\
    } while (0)

/*
 * Locate area index for *incmp
 * Set area lower / upper bound 0 if *incmp < min(b), (BINDEX_K - 1) if *incmp > max(b)
 * Set position lower / upper bound of the absolute position for *incmp
*/
void
BDXlocate(BAT* b, void *incmp, int *al, int *au, BUN *pl, BUN *pu)
{
	Bindex *bindex = b->tbindex;
	area_t* area = (area_t*)(bindex->area);
	const BUN * restrict rv = (BUN *)(bindex->rv);

	switch (ATOMbasetype(b->ttype)) {
	case TYPE_bte:
		LOCATE_POS(bte);
		break;
	case TYPE_sht:
		LOCATE_POS(sht);
		break;
	case TYPE_int:
		LOCATE_POS(int);
		break;
	case TYPE_lng:
		LOCATE_POS(lng);
		break;
#ifdef HAVE_HGE
	case TYPE_hge:
		LOCATE_POS(hge);
		break;
#endif
	default:
		assert(0);
		break;
	}
}

/*
 * If anti is set, copy the bitwise-not fv to sv
 */
void
BDXcopyfv(BAT *b, int fvi, bool anti, BUN bmn)
{
	bitmap_t *sv = (bitmap_t *)b->tbindex->sv;
	bitmap_t *fv = (bitmap_t *)b->tbindex->fv + fvi * bmn;

	assert(fvi >= -1);
	assert(fvi < BINDEX_K);

	if ((fvi == -1 && !anti) || ((fvi == BINDEX_K - 1) && anti)) {
		memset(sv, 0, bmn * sizeof(bitmap_t));
	} else if ((fvi == -1 && anti) || ((fvi == BINDEX_K - 1) && !anti)) {
		memset(sv, 0xff, bmn * sizeof(bitmap_t));
	} else {
		for (BUN i = 0; i < bmn; i++)
			sv[i] = anti ? (~fv[i]) : fv[i];
	}
	sv[bmn - 1] &= ~0UL >> (bmn * BITMAPWIDTH - BATcount(b));
}

/*
 * If anti is set, copy the bitwise-not fv.
 * If andop is set, merge with bitwise-and, otherwise merge with bitwise-or
 */
void
BDXmergefv(BAT *b, int fvi1, int fvi2, bool anti1, bool anti2, bool andop, BUN bmn)
{
	bitmap_t *sv = (bitmap_t *)b->tbindex->sv;
	bitmap_t *fv1 = (bitmap_t *)b->tbindex->fv + fvi1 * bmn;
	bitmap_t *fv2 = (bitmap_t *)b->tbindex->fv + fvi2 * bmn;
	assert(fvi1 >= -1);
	assert(fvi1 < BINDEX_K);
	assert(fvi2 >= -1);
	assert(fvi2 < BINDEX_K);

	bool zero1 = (fvi1 == -1 && !anti1) || ((fvi1 == BINDEX_K - 1) && anti1);
	bool zero2 = (fvi2 == -1 && !anti2) || ((fvi2 == BINDEX_K - 1) && anti2);
	bool one1 = (fvi1 == -1 && anti1) || ((fvi1 == BINDEX_K - 1) && !anti1);
	bool one2 = (fvi2 == -1 && anti2) || ((fvi2 == BINDEX_K - 1) && !anti2);
	assert(!(zero1 && one1));
	assert(!(zero2 && one2));

	if ((andop && (zero1 || zero2)) || (!andop && (zero1 && zero2))) {
		/* all 0 */
		memset(sv, 0, bmn * sizeof(bitmap_t));
	} else if ((andop && (one1 && one2)) || (!andop && (one1 || one2))) {
		/* all 1 */
		memset(sv, 0xff, bmn * sizeof(bitmap_t));
	} else if ((andop && one1) || (!andop && zero1)) {
		/* copy fv2 */
		for (BUN i = 0; i < bmn; i++)
			sv[i] = anti2 ? (~fv2[i]) : fv2[i];
	} else if ((andop && one2) || (!andop && zero2)) {
		/* copy fv1 */
		for (BUN i = 0; i < bmn; i++)
			sv[i] = anti1 ? (~fv1[i]) : fv1[i];
	} else {
		/* merge fv1 and fv2 */
		for (BUN i = 0; i < bmn; i++) {
			bitmap_t bitset1 = anti1 ? (~fv1[i]) : fv1[i];
			bitmap_t bitset2 = anti2 ? (~fv2[i]) : fv2[i];
			sv[i] = andop ? (bitset1 & bitset2) : (bitset1 | bitset2);
		}
	}
	sv[bmn - 1] &= ~0UL >> (bmn * BITMAPWIDTH - BATcount(b));
}

/*
 * Implementation of BDXrefinesv
 */
static void
BDXrefinesv_(void *sv, BUN *rv, BUN s, BUN e)
{
	bitmap_t *sv_ = (bitmap_t *)sv;
	BUN *cur = rv + s;

	while (cur < rv + e) {
		BUN pos = *cur++;
		sv_[pos / BITMAPWIDTH] ^= 1UL << (pos % BITMAPWIDTH);
	}
}

/*
 * Refine scratch vector fv according to rv[s, e)
 */
void
BDXrefinesv(BAT *b, BUN s, BUN e)
{
	return BDXrefinesv_(b->tbindex->sv, b->tbindex->rv, s, e);
}

/*
 * Merge 2 BAT's bitmap
 */
BAT *
BDXmergesv(BAT *b, BAT *s)
{
	assert(b->tbindex);
	BAT *bn = BATdense(0, 0, 0);
	BATsetcount(bn, b->batCount);
	bn->tbm = b->tbindex->sv;
	bn->tbmn = b->tbindex->bmn;

	// a complete column
	if (!s || (s && !s->tbm))
		return bn;

	// bitwise-and on 2 bitmaps
	assert(s->tbmn == bn->tbmn);
	BUN bmn = bn->tbmn;
	bitmap_t *bnbm = (bitmap_t *)bn->tbm;
	bitmap_t *sbm = (bitmap_t *)s->tbm;

	for (BUN i = 0; i < bmn; i++)
		bnbm[i] &= sbm[i];

	return bn;
}

/*
 * Construct a BAT * from it's sv
 */
BAT *
BDXsv2bat(BAT *b)
{
	void *sv = b->tbindex->sv;
	BUN bmn = b->tbindex->bmn;
	BAT *bn;

	bn = COLnew(0, TYPE_oid, bmn * BITMAPWIDTH, TRANSIENT);
	if (bn == NULL)
		return NULL;

	oid *rbn = (oid *) Tloc((bn), 0);
	BUN cnt = 0;

	for (BUN i = 0; i < bmn; i++) {
		bitmap_t bitset = *((bitmap_t *)sv + i);
		while(bitset) {
			int ntz = __builtin_ctzll(bitset);
			rbn[cnt++] = i * BITMAPWIDTH + ntz;
			bitset ^= (1UL << ntz);
		}
	}

	assert(bn->batCapacity >= cnt);
	BATsetcount(bn, cnt);

	bn->tsorted = true;
	bn->trevsorted = bn->batCount <= 1;
	bn->tkey = true;
	bn->tseqbase = cnt == 0 ? 0 : cnt == 1 || cnt == b->batCount ? b->hseqbase : oid_nil;
	bn->tnil = false;
	bn->tnonil = true;

	return bn;
}

#define CREATE_BINDEX(TYPE)	\
	do {	\
		BUN cnt = BATcount(b);	\
		TYPE *restrict scol = (TYPE *) Tloc(s, 0);	\
		BUN *restrict ocol = (BUN *) Tloc(o, 0);	\
		BUN arean = cnt / BINDEX_K;	\
		/* build area */	\
		for (BUN i = 0; i < BINDEX_K; i++) {	\
			BUN p = i * arean;	\
			area[i].val.TYPE##val = scol[p];	\
			area[i].vtype = tpe;	\
			area[i].pos = p;	\
			area[i].len = (i == BINDEX_K - 1) ? (cnt - p) : arean;	\
		}	\
		/* build refine vector */	\
		if (BATtdense(o)) {	\
			/* the input is sorted */	\
			for (BUN i = 0; i < cnt; i++)	\
				rv[i] = i;	\
		} else {	\
			assert(ocol);	\
			memcpy(rv, ocol, cnt * sizeof(BUN));	\
		}	\
		/* build filter vector */	\
		for (BUN i = 0; i < BINDEX_K - 1; i++) {	\
			bitmap_t * fvi = (bitmap_t *)fv + i * bmn;	\
			bitmap_t * fvp = (bitmap_t *)fv + (i - 1) * bmn;	\
			if (i == 0)	\
				memset(fvi, 0, bmn * sizeof(bitmap_t));	\
			else	\
				memcpy(fvi, fvp, bmn * sizeof(bitmap_t));	\
			BDXrefinesv_(fvi, rv, area[i].pos, area[i].pos + area[i].len);	\
		}	\
		/* clear scratch vector */	\
		memset(sv, 0, bmn * sizeof(bitmap_t));	\
	} while (0);

static void
bindex_create(BAT* b, BAT* s, BAT* o, void *fv, void *sv, area_t *area, BUN *rv, BUN bmn)
{
	int tpe = ATOMbasetype(b->ttype);
	switch (tpe) {
	case TYPE_bte:
		CREATE_BINDEX(bte);
		break;
	case TYPE_sht:
		CREATE_BINDEX(sht);
		break;
	case TYPE_int:
		CREATE_BINDEX(int);
		break;
	case TYPE_lng:
		CREATE_BINDEX(lng);
		break;
#ifdef HAVE_HGE
	case TYPE_hge:
		CREATE_BINDEX(hge);
		break;
#endif
	default:
		/* should never reach here */
		assert(0);
	}
}

/* Check whether we have bindex on b (and return true if we do).  It
 * may be that the bindex were made persistent, but we hadn't seen
 * that yet, so check the file system.  This also returns true if b is
 * a view and there are bindex on b's parent.
 *
 * Note that the b->tbindex pointer can be NULL, meaning there are
 * no bindex; (Bindex *) 1, meaning there are no bindex loaded,
 * but they may exist on disk; or a valid pointer to loaded bindex.
 * These values are maintained here, in the BDXdestroy and BDXfree
 * functions, and in BBPdiskscan during initialization. */
// TODO: bindex BBPdiskscan
bool
BATcheckbindex(BAT *b)
{
	bool ret;

	if (VIEWtparent(b)) {
		assert(b->tbindex == NULL);
		b = BBPdescriptor(VIEWtparent(b));
	}

	if (b->tbindex == (Bindex *) 1) {
		MT_lock_set(&GDKbindexLock(b->batCacheid));
		if (b->tbindex == (Bindex *) 1) {
			Bindex *bindex;
			const char *nme = BBP_physical(b->batCacheid);

			b->tbindex = NULL;
			if ((bindex = GDKzalloc(sizeof(Bindex))) != NULL &&
			    (bindex->bindex.farmid = BBPselectfarm(b->batRole, b->ttype, bindexheap)) >= 0) {
				int fd;

				stpconcat(bindex->bindex.filename, nme, ".bindex", NULL);
				/* check whether a persisted bindex index
				 * can be found */
				if ((fd = GDKfdlocate(bindex->bindex.farmid, nme, "rb", "tbindex")) >= 0) {
					size_t hdata[BINDEX_HEADER_SIZE];
					struct stat st;

					if (read(fd, hdata, sizeof(hdata)) == sizeof(hdata) &&
						hdata[0] & ((size_t) 1 << 16) &&
						((hdata[0] & 0xFF00) >> 8) == BINDEX_VERSION &&
						hdata[2] == (size_t) BATcount(b) &&
						fstat(fd, &st) == 0 &&
						st.st_size >= (off_t) (bindex->bindex.size =
									bindex->bindex.free =
									BINDEX_HEADER_SIZE * SIZEOF_SIZE_T + /* extra info */
									(BINDEX_K - 1) * hdata[1] * sizeof(bitmap_t) + /* filter vector[K - 1] */
									1 * hdata[1] * sizeof(bitmap_t) + /* scratch vector */
									BINDEX_K * sizeof(area_t) + /* area[K] */
									hdata[2] * sizeof(BUN) /* refine vector */
									) &&
						HEAPload(&bindex->bindex, nme, "tbindex", false) == GDK_SUCCEED) {
						/* usable */
						bindex->bmn = (BUN) hdata[1];
						bindex->fv = bindex->bindex.base + BINDEX_HEADER_SIZE * SIZEOF_SIZE_T;
						bindex->area = (void *) ((bitmap_t *)bindex->fv + (BINDEX_K - 1) * bindex->bmn);
						bindex->rv = (void *) ((area_t *)bindex->area + BINDEX_K);
						close(fd);
						bindex->bindex.parentid = b->batCacheid;
						b->tbindex = bindex;
						ALGODEBUG fprintf(stderr, "#BATcheckbindex(" ALGOBATFMT "): reusing persisted bindex\n", ALGOBATPAR(b));
						MT_lock_unset(&GDKbindexLock(b->batCacheid));

						return true;
					}
					close(fd);
					/* unlink unusable file */
					GDKunlink(bindex->bindex.farmid, BATDIR, nme, "tbindex");
				}
			}
			GDKfree(bindex);
			GDKclrerr();	/* we're not currently interested in errors */
		}
		MT_lock_unset(&GDKbindexLock(b->batCacheid));
	}
	ret = b->tbindex != NULL;
	ACCELDEBUG if (ret) fprintf(stderr, "#BATcheckbindex(" ALGOBATFMT "): already has bindex\n", ALGOBATPAR(b));
	return ret;
}

static void
BATbdxsync(void *arg)
{
	BAT *b = arg;
	Bindex *bindex;
	int fd;
	lng t0 = 0;
	const char *failed = " failed";

	ACCELDEBUG t0 = GDKusec();

	MT_lock_set(&GDKbindexLock(b->batCacheid));
	if ((bindex = b->tbindex) != NULL) {
		Heap *hp = &bindex->bindex;
		if (HEAPsave(hp, hp->filename, NULL) == GDK_SUCCEED) {
			if (hp->storage == STORE_MEM) {
				if ((fd = GDKfdlocate(hp->farmid, hp->filename, "rb+", NULL)) >= 0) {
					/* add version number */
					((size_t *) hp->base)[0] |= (size_t) BINDEX_VERSION << 8;
					/* sync-on-disk checked bit */
					((size_t *) hp->base)[0] |= (size_t) 1 << 16;
					if (write(fd, hp->base, SIZEOF_SIZE_T) >= 0) {
						failed = ""; /* not failed */
						if (!(GDKdebug & NOSYNCMASK)) {
#if defined(NATIVE_WIN32)
							_commit(fd);
#elif defined(HAVE_FDATASYNC)
							fdatasync(fd);
#elif defined(HAVE_FSYNC)
							fsync(fd);
#endif
						}
						hp->dirty = false;
					} else {
						failed = " write failed";
						perror("write hash");
					}
					close(fd);
				}
			} else {
				/* add version number */
				((size_t *) hp->base)[0] |= (size_t) BINDEX_VERSION << 8;
				/* sync-on-disk checked bit */
				((size_t *) hp->base)[0] |= (size_t) 1 << 16;
				if (!(GDKdebug & NOSYNCMASK) &&
					MT_msync(hp->base, SIZEOF_SIZE_T) < 0) {
					failed = " sync failed";
					((size_t *) hp->base)[0] &= ~((size_t) BINDEX_VERSION << 8);
				} else {
					hp->dirty = false;
					failed = ""; /* not failed */
				}
			}
			ACCELDEBUG fprintf(stderr, "#BATbdxsync(" ALGOBATFMT "): "
					  "bindex persisted "
					  "(" LLFMT " usec)%s\n", ALGOBATPAR(b),
					  GDKusec() - t0, failed);
		}
	}
	MT_lock_unset(&GDKbindexLock(b->batCacheid));
	BBPunfix(b->batCacheid);
}

gdk_return
BATbindex(BAT *b)
{
	BAT *s1 = NULL, *s2 = NULL;
	BAT *o = NULL;
	Bindex *bindex;
	lng t0 = 0;
	BUN cnt;

	/* we only create bindex for types that look like types we know */
	switch (ATOMbasetype(b->ttype)) {
	case TYPE_bte:
	case TYPE_sht:
	case TYPE_int:
	case TYPE_lng:
#ifdef HAVE_HGE
	case TYPE_hge:
#endif
	break;
	default:		/* type not supported */
		/* doesn't look enough like base type: do nothing */
		GDKerror("BATbindex: unsupported type\n");
		return GDK_FAIL;
	}

	BATcheck(b, "BATbindex", GDK_FAIL);

	/* not support nill value */
	if (b->tnil)
		return GDK_FAIL;

	/* not enough records */
	cnt = BATcount(b);
	if (cnt < BINDEX_K)
		return GDK_FAIL;

	if (BATcheckbindex(b))
		return GDK_SUCCEED;

	if (VIEWtparent(b)) {
		/* views always keep null pointer and need to obtain
		 * the latest bindex from the parent at query time */
		s2 = b;		/* remember for ACCELDEBUG print */
		b = BBPdescriptor(VIEWtparent(b));
		assert(b);
		if (BATcheckbindex(b))
			return GDK_SUCCEED;
	}
	MT_lock_set(&GDKbindexLock(b->batCacheid));
	ACCELDEBUG t0 = GDKusec();
	if (b->tbindex == NULL) {
		const char *nme = BBP_physical(b->batCacheid);

		ACCELDEBUG {
			if (s2)
				fprintf(stderr, "#BATbindex(b=" ALGOBATFMT
					"): creating bindex on parent "
					ALGOBATFMT "\n",
					ALGOBATPAR(s2), ALGOBATPAR(b));
			else
				fprintf(stderr, "#BATbindex(b=" ALGOBATFMT
					"): creating bindex\n",
					ALGOBATPAR(b));
		}
		s2 = NULL;

		bindex = GDKzalloc(sizeof(Bindex));
		if (bindex == NULL) {
			MT_lock_unset(&GDKbindexLock(b->batCacheid));
			return GDK_FAIL;
		}
		stpconcat(bindex->bindex.filename, nme, ".tbindex", NULL);
		bindex->bindex.farmid = BBPselectfarm(b->batRole, b->ttype,
							   bindexheap);

		if (BATsort(&s1, &o, NULL, b, NULL, NULL, false, false, false) != GDK_SUCCEED) {
			MT_lock_unset(&GDKbindexLock(b->batCacheid));
			GDKfree(bindex);
			return GDK_FAIL;
		}

		assert(BATcount(b) == BATcount(o));
		BUN bmn = ((cnt + BITMAPWIDTH - 1) / BITMAPWIDTH);
		bindex->bmn = bmn;

		/* The heap we create here consists of five parts:
		 * fv, filter vectors[K - 1];
		 * sv, scratch vector;
		 * area, area[K];
		 * rv, refine vector;
		 * bmn, bitmap length;
		 * In addition, we add some housekeeping entries at
		 * the start so that we can determine whether we can
		 * trust the bindex when encountered on startup (including
		 * a version number -- CURRENT VERSION is 1). */
		if (HEAPalloc(&bindex->bindex,
			    BINDEX_HEADER_SIZE * SIZEOF_SIZE_T + /* extra info */
				(BINDEX_K - 1) * bmn * sizeof(bitmap_t) + /* filter vector[K - 1] */
				1 * bmn * sizeof(bitmap_t) + /* scratch vector */
				BINDEX_K * sizeof(area_t) + /* area[K] */
				cnt * sizeof(BUN), /* refine vector */
				1) != GDK_SUCCEED) {
			MT_lock_unset(&GDKbindexLock(b->batCacheid));
			GDKfree(bindex);
			GDKerror("#BATbindex: memory allocation error");
			BBPunfix(s1->batCacheid);
			BBPunfix(o->batCacheid);
			return GDK_FAIL;
		}
		bindex->fv = bindex->bindex.base + BINDEX_HEADER_SIZE * SIZEOF_SIZE_T;
		bindex->sv = (void *) ((bitmap_t *)bindex->fv + (BINDEX_K - 1) * bmn);
		bindex->area = (void *) ((bitmap_t *)bindex->sv + bmn);
		bindex->rv = (void *) ((area_t *)bindex->area + BINDEX_K);

		bindex_create(b, s1, o, bindex->fv, bindex->sv, (area_t *)bindex->area, bindex->rv, bindex->bmn);

		bindex->bindex.free = (size_t) ((char *) ((BUN *) bindex->rv + cnt) - bindex->bindex.base);
		/* add info to heap for when they become persistent */
		((size_t *) bindex->bindex.base)[0] = (size_t)0;
		((size_t *) bindex->bindex.base)[1] = (size_t)bindex->bmn;
		((size_t *) bindex->bindex.base)[2] = (size_t)BATcount(b);
		bindex->bindex.parentid = b->batCacheid;
		b->tbindex = bindex;
		if (BBP_status(b->batCacheid) & BBPEXISTING &&
		    !b->theap.dirty) {
			MT_Id tid;
			BBPfix(b->batCacheid);
			char name[16];
			snprintf(name, sizeof(name), "bdxsync%d", b->batCacheid);
			if (MT_create_thread(&tid, BATbdxsync, b,
					     MT_THR_DETACHED, name) < 0)
				BBPunfix(b->batCacheid);
		}
	}

	ACCELDEBUG fprintf(stderr, "#BATbindex(%s): bindex construction " LLFMT " usec\n", BATgetId(b), GDKusec() - t0);
	MT_lock_unset(&GDKbindexLock(b->batCacheid));

	// TODO: bindex dead lock here?

	return GDK_SUCCEED;
}

lng
BDXbindexsize(BAT *b)
{
	lng sz = 0;
	if (b->tbindex && b->tbindex != (Bindex *) 1) {
		sz += (BINDEX_K - 1) * b->tbindex->bmn * sizeof(bitmap_t);
		sz += b->tbindex->bmn * sizeof(bitmap_t);
		sz += BINDEX_K * sizeof(area_t);
		sz += BATcount(b) * sizeof(BUN);
	}
	return sz;
}

static void
BDXremove(BAT *b)
{
	Bindex *bindex;

	assert(b->tbindex != NULL);
	assert(!VIEWtparent(b));

	if ((bindex = b->tbindex) != NULL) {
		b->tbindex = NULL;

		if ((GDKdebug & ALGOMASK) &&
		    * (size_t *) bindex->bindex.base & (1 << 16))
			fprintf(stderr, "#BDXremove: removing persisted bindex\n");
		if (HEAPdelete(&bindex->bindex, BBP_physical(b->batCacheid),
			       "tbindex") != GDK_SUCCEED)
			IODEBUG fprintf(stderr, "#BDXremove(%s): bindex heap\n", BATgetId(b));

		GDKfree(bindex);
	}
}

void
BDXdestroy(BAT *b)
{
	if (b && b->tbindex) {
		MT_lock_set(&GDKbindexLock(b->batCacheid));
		if (b->tbindex == (Bindex *) 1) {
			b->tbindex = NULL;
			GDKunlink(BBPselectfarm(b->batRole, b->ttype, bindexheap),
				  BATDIR,
				  BBP_physical(b->batCacheid),
				  "tbindex");
		} else if (b->tbindex != NULL && !VIEWtparent(b))
			BDXremove(b);
		MT_lock_unset(&GDKbindexLock(b->batCacheid));
	}
}

/* free the memory associated with the bindex, do not remove the
 * heap files; indicate that bindex are available on disk by setting
 * the bindex pointer to 1 */
void
BDXfree(BAT *b)
{
	Bindex *bindex;

	if (b && b->tbindex) {
		assert(b->batCacheid > 0);
		MT_lock_set(&GDKbindexLock(b->batCacheid));
		bindex = b->tbindex;
		if (bindex != NULL && bindex != (Bindex *) 1) {
			b->tbindex = (Bindex *) 1;
			if (!VIEWtparent(b)) {
				HEAPfree(&bindex->bindex, false);
				GDKfree(bindex);
			}
		}
		MT_lock_unset(&GDKbindexLock(b->batCacheid));
	}
}

#ifndef NDEBUG
/* never called, useful for debugging */

static void
printfv(void *fvi_, BUN bmn)
{
#define PRTBITMAP(bitmap)	\
	do {	\
		for (int i = BITMAPWIDTH - 1; i >= 0 ; i--)	\
			fprintf(stderr, "%d", !!(bitmap & (1UL << i)));	\
	} while (0);

	bitmap_t * fvi = (bitmap_t *)fvi_;
	for (BUN j = 0; j < bmn; j++) {
		PRTBITMAP(fvi[j]);
		fprintf(stderr, " | ");
	}
	fprintf(stderr, "\n");
}

static void
printarea(area_t *area, int tpe)
{
#define PRTVAL(FMT, TYPE) fprintf(stderr, "area->val = " FMT "\n", area->val.TYPE);

	fprintf(stderr, "area->pos = " BUNFMT "\t", area->pos);
	fprintf(stderr, "area->len = " BUNFMT "\t", area->len);
	switch (tpe) {
	case TYPE_bte:
		PRTVAL("%d", bteval);
		break;
	case TYPE_sht:
		PRTVAL("%d", shtval);
		break;
	case TYPE_int:
		PRTVAL("%d", intval);
		break;
	case TYPE_lng:
		PRTVAL("%ld", lngval);
		break;
#ifdef HAVE_HGE
	case TYPE_hge:
		// TODO: bindex display hge value
		// PRTVAL(, hgeval);
		break;
#endif
	default:
		/* should never reach here */
		assert(0);
	}
	fprintf(stderr, "\n");
}

static void
printrv(BUN *rv, BUN cnt)
{
	for (BUN i = 0; i < cnt; i++) {
		fprintf(stderr, BUNFMT " ", rv[i]);
		if ((i + 1) % 16 == 0)
			fprintf(stderr, "\n");
	}
	fprintf(stderr,"\n");
}

void BDXprint(BAT *b)
{
	Bindex *bindex;
	BUN cnt, bmn;
	bitmap_t *fv;
	bitmap_t *sv;
	area_t *area;
	BUN *rv;

	if (!BATcheckbindex(b)) {
		fprintf(stderr, "no bindex\n");
		return;
	}

	bindex = b->tbindex;
	cnt = BATcount(b);
	bmn = bindex->bmn;
	fv = (bitmap_t *)bindex->fv;
	sv = (bitmap_t *)bindex->sv;
	area = (area_t *)bindex->area;
	rv = bindex->rv;

	fprintf(stderr, "cnt = " BUNFMT ", bmn = " BUNFMT "\n", cnt, bmn);

	fprintf(stderr, "fv\n");
	for (int i = 0; i < BINDEX_K - 1; i++) {
		fprintf(stderr, "fv[%d]\n", i);
		bitmap_t *fvi = fv + i * bmn;
		printfv(fvi, bmn);
	}
	fprintf(stderr,"\n");

	fprintf(stderr, "area\n");
	for (int i = 0; i < BINDEX_K; i++) {
		fprintf(stderr, "area[%d]\n", i);
		printarea(&area[i], ATOMbasetype(b->ttype));
	}
	fprintf(stderr,"\n");

	fprintf(stderr, "rv\n");
	printrv(rv, cnt);
	fprintf(stderr,"\n");

	fprintf(stderr, "sv\n");
	printfv(sv, bmn);
	fprintf(stderr,"\n");
}

#endif
