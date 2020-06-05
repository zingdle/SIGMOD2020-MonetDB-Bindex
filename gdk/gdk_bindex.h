/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * Jiading Guo (zingdle@gmail.com)
 */

#ifndef GDK_BINDEX_H
#define GDK_BINDEX_H

typedef uint64_t bitmap_t;

#define BITMAPWIDTH (sizeof(bitmap_t) << 3)

#define BINDEX_K 128

/*
 * the area struct
 */
typedef struct {
    union {			/* storage is first in the record */
		bte bteval;
		sht shtval;
		int intval;
		lng lngval;
#ifdef HAVE_HGE
		hge hgeval;
#endif
	} val;
	int vtype;
    BUN pos;
    BUN len;
} area_t;

typedef enum {
	op_le = 0,	/* v <= vh */
	op_ge,		/* v >= vl */
	op_bt,		/* v >= vl && v <= vh */
	op_ab,		/* (v <= vl || v >= vh) */
} op_t;

#endif  /* GDK_BINDEX_H */
