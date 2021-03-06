//
//  Copyright (c) 2014, 2015, 2016 Beijing Genomics Institution (BGI)
//  shiquan@genomics.cn
//

#include "commons.h"
#include "bedutil.h"
#include "ksort.h"
#include "kseq.h"
#include "kstring.h"
#include "number.h"

KSORT_INIT_GENERIC(uint64_t)

KSTREAM_INIT(gzFile, gzread, 8192)

#define LIDX_SHIFT 13
#define swapvalue(a, b, t) {t c = a; a = b; b = c;}

static void reg_sort(reglist_t *reg)
{
#ifdef _DEBUG_MODE
    debug("[%s]", __func__);
#endif
    
    if (reg->flag & BD_IS_SORTED) {
	debug("Already sorted!");
	return;
    }
    ks_introsort(uint64_t, reg->m, reg->a);
    reg->flag |= BD_IS_SORTED;
}

static int get_id(bedaux_t *reg, char *s)
{
    int i;
    for ( i = 0; i < reg->n_seq; ++i )
    {
	if ( ! strcmp(s, reg->seq_names[i]) ) return i;
    }
    return -1;
}

/*
 * read bed file (0-based chromosome begin end) and tsv file (1-based chromosome posotion)
 */
static void bed_read(const char *fn, bedaux_t * reg, int right_flank, int left_flank, int *is_error, int trim_tag) 
{
#ifdef _DEBUG_MODE
    debug("[%s] %s", __func__,fn);
#endif

    assert ( right_flank >= 0 && left_flank >=0 );
    assert ( reg->alloced >= reg->n_files );
    if ( reg->alloced == reg->n_files ) {
	reg->alloced = reg->alloced == 0 ? 2 : reg->alloced + 10;
	reg->filenames = (char**)realloc(reg->filenames, reg->alloced*sizeof(char*));
	reg->hfiles = (bedfile_t**)realloc(reg->hfiles, reg->alloced*sizeof(bedfile_t*));
    }
    int nfile = reg->n_files;
    int i;
    for (i = 0; nfile > 0 && i < nfile; i++) {
	if ( reg->filenames[i] && !strcmp(fn, reg->filenames[i]) ) {
	    warnings("File %s load more than once. Skip...", fn);
	    return;
	}
    }
    gzFile fp = safe_gzopen(fn); // maybe we may skip these NULL files instead of errabort
    reg->filenames[nfile] = strdup(fn);
    reg->hfiles[nfile] = (bedfile_t*)malloc(sizeof(bedfile_t));
    reg->hfiles[nfile]->is_empty = BD_IS_EMPTY;
    reg->hfiles[nfile]->reg = kh_init(reg);
    reg->n_files++;
    kstream_t * ks;
    int dret;
    uint32_t length = 0;
    kstring_t *str;
    regHash_t *reghash = reg->hfiles[nfile]->reg;
    str = (kstring_t*)needmem(sizeof(kstring_t));
    ks = ks_init(fp);
    int line = 0;
    while (ks_getuntil(ks, 2, str, &dret) >= 0) {
	int32_t beg = -1, end = -1;
	line++;
	if (isNull(str->s) || str->s[0] == '\n') {
	    //if (dret != '\n') while ((ks_getc(ks)) > 0 && dret != '\n');
	    warnings("%s: line %d is empty! skip... ", fn, line);
	    continue;
	}

	if (str->s[0] == '#') {
	    //while (dret != '\n' && ks_getuntil(ks, 0, str, &dret) > 0);
	    continue;
	}
	int nfields = 0;
	int *splits = ksplit(str, 0, &nfields);
	if (splits == NULL)
	    continue;

	if (nfields < 2) {
	    free(splits);
	    continue;
	}
	char *name = str->s + splits[0];
	char *tmp = str->s + splits[1];
	if ( check_num_likely(tmp) ) {
	    beg = str2int(tmp);
        } else {
	    warnings("%s: line %d is malformed! skip... ", fn, line);
	    continue;
	}
	if (nfields > 2 && check_num_likely(str->s+splits[2])) {
	    end = str2int(str->s + splits[2]);
	}
	free(splits);
	khiter_t k;
	k= kh_get(reg, reghash, name);
	if (k == kh_end(reghash)) {
	    int ret;
	    reglist_t *b;
	    b = (reglist_t*)needmem(sizeof(reglist_t));
	    b->sorted = BD_IS_UNSORT;
	    b->id = get_id(reg, name);
	    if (b->id == -1 ) {
		b->id = reg->n_seq;
		reg->n_seq++;
		reg->seq_names = (char**)realloc(reg->seq_names, reg->n_seq*sizeof(char*));
		reg->seq_names[b->id] = strdup(name);
	    }
	    k = kh_put(reg, reghash, reg->seq_names[b->id], &ret);
	    kh_val(reghash, k) = b;      
	}
	reglist_t *p = kh_val(reghash, k);
	if (beg == end) *is_error = 1;
	if (end == -1) { // treat as tsv file 
	    // this is different from bed_read in bedidx.c
	    end = beg;
	    beg = beg < 1 ? 0 : beg-1;
	}
	if (end < beg) swapvalue(beg, end, uint32_t);
	length += end - beg;
	if ( trim_tag == BD_IS_TRIM) {
	    if ( end <= right_flank + left_flank + beg ) {
		warnings("region %d-%d is too short to trim! Skip...", beg, end);
	    } else {
		beg += left_flank;
		end -= right_flank;
	    }
	} else {
	    beg -= left_flank;
	    end += right_flank;
	}
	if (beg < 0) beg = 0; // 0baesd, first position
	if (p->m == p->n) {
	    p->n = p->n ? p->n * 2 : 4;
	    p->a = (uint64_t*)realloc(p->a, p->n*sizeof(uint64_t));
	}
	p->a[p->m++] = (uint64_t)beg<<32 | (uint32_t)end;
	str->l = 0;
    }

    bedfile_t *bed = reg->hfiles[nfile];
    bed->region = 0;
    bed->length = 0;
    khiter_t k;
    uint32_t regions = 0;
    for ( k = 0; k < kh_end(bed->reg); ++k ) {
	if (kh_exist(bed->reg, k)) {
	    regions += kh_val(bed->reg, k)->m;
	}
    }
    if ( regions ) {
	reg->hfiles[nfile]->is_empty = BD_IS_HAVE;
	bed->region = regions;
    } else {
	reg->hfiles[nfile]->is_empty = BD_IS_EMPTY;
    }
    bed->length += length;
    reg->length += bed->length;
    reg->region += bed->region;
    //reg->region += 
    ks_destroy(ks);
    gzclose(fp);
    freemem(str->s); freemem(str);
}

static void clear_reg(reglist_t *reg, bedvoid_destroy func)
{
#ifdef _DEBUG_MODE
    debug("[%s]", __func__);
#endif

    if ( reg == NULL ) return;
    if ( reg->n == 0 ) return;
    freemem(reg->a);
    if (reg->count) freemem(reg->count);
    freemem(reg->idx);
    func(reg->data);
}

static void bed_clear(bedaux_t * reg, bedvoid_destroy func) 
{
#ifdef _DEBUG_MODE
    debug("[%s]", __func__);
#endif

    int i, j;
    khiter_t k;
    for ( i = 0; i < reg->n_seq; ++i ) 
    {
	for ( j = 0; j < reg->n_files; ++j)
	{
	    if ( reg->hfiles[j] == NULL || reg->hfiles[j]->is_empty == BD_IS_EMPTY ) continue;
	    regHash_t *rgh = reg->hfiles[j]->reg;
	    k = kh_get(reg, rgh, reg->seq_names[i]);
	    if ( k != kh_end(rgh) ) {
		reglist_t *bed = kh_val(rgh, k);
		clear_reg(bed, func);
		freemem(bed);
		kh_del(reg, rgh, k);
	    }
	}
	freemem(reg->seq_names[i]);
    }
    for ( i = 0; i < reg->n_files; ++i)
    {
	regHash_t *rgh = reg->hfiles[i]->reg;
	kh_destroy(reg, rgh);
	freemem(reg->filenames[i]);
	freemem(reg->hfiles[i]);
    }
  
    freemem(reg->filenames);
    freemem(reg->hfiles);
    freemem(reg->seq_names);
}
static void bed_destroy(bedaux_t * reg, bedvoid_destroy func)
{
#ifdef _DEBUG_MODE
    debug("[%s]", __func__);
#endif

    bed_clear(reg, func);
    freemem(reg);
}

void destroy_reg(void *data)
{
#ifdef _DEBUG_MODE
    debug("[%s]", __func__);
#endif
    
    reglist_t *bed = (reglist_t*)data;
    if (bed==NULL) return;
    clear_reg(bed, destroy_void);
    freemem(bed);
}

void destroy_void(void *data)
{
    return;
}

static reglist_t * reg_add(reglist_t **beds, int n_beds)
{
#ifdef _DEBUG_MODE
    debug("[%s] n_beds : %d", __func__, n_beds);
#endif
    if (n_beds == 0)
	return NULL;
    assert(beds != NULL);
    int i;
    //for(i=0; i<n_beds && beds[i] == NULL; ++i);
    //if ( i == n_beds-1 ) return NULL;
    
    reglist_t * bed = (reglist_t*)calloc(1,sizeof(reglist_t));
    bed->sorted = BD_IS_UNSORT;
    bed->id = -1;
    bed->a = NULL;
    bed->n = bed->m = 0;
    bed->count = NULL;
    for (i=0; i < n_beds; ++i )
    {	
	if ( beds[i] == NULL || beds[i]->m == 0 ) continue;
	if (bed->id == -1)
	    bed->id = beds[i]->id;
	else if ( beds[i]->id != bed->id )
	    errabort("[%s]: beds[%d] != bed->id : %d\t%d", __func__, i, beds[i]->id, bed->id);

	if ( bed->n <= bed->m + beds[i]->m ) {
	    bed->n = bed->m + beds[i]->m;
	    assert(bed->n >0);
	    kroundup32(bed->n);
	    bed->a = (uint64_t*)realloc(bed->a, bed->n *sizeof(uint64_t));
	}
	memcpy(bed->a + bed->m, beds[i]->a, beds[i]->m * sizeof(uint64_t));
	bed->m += beds[i]->m;
	bed->l_reg += beds[i]->l_reg;
    }
    if (bed->m == 0) {
	if (bed->n)
	    free(bed->a);
	free(bed);
	return NULL;
    }
    bed->count = (uint32_t*)malloc((bed->m+1)*sizeof(uint32_t));
    return bed;
}

static int check_reg_by_chromosome_length(reglist_t * bed)
{
    if ( bed == NULL ) return -1;
    if ( bed->sorted == BD_IS_UNSORT ) return -2;
    if ( bed->l_chr && bed->l_chr > maxlen(bed) ) return 0;
    int i;
    for ( i = bed->m; i >= 0; i-- )
    {
	if ( bed->a[i]>>32 >= bed->l_chr )
	    continue;
	else if ( (uint32_t)bed->a[i] >= bed->l_chr ) {
	    bed->a[i] = (bed->a[i] & (uint32_t) 0) | bed->l_chr;
	} else
	    errabort("[%s] Unsorted bed!", __func__);
    }
    return 1;
}

static void regcore_merge(reglist_t *bed) 
{
#ifdef _DEBUG_MODE
    debug("[%s]", __func__);
#endif

    assert(bed != NULL);

    if ( bed->sorted & BD_IS_MERGED ) {
	errabort("bed is already merged!");
	return;
    }
    uint32_t length = 0;
    if (bed->m == 0) {
	errabort("bed is empty m==0");
	return;
    }
    int i, m = 0;
    if (bed->count == NULL) 
	bed->count = (uint32_t*)calloc(bed->m, sizeof(uint32_t));

    for ( i = 0; i < bed->m; ++i)
	bed->count[i] = 1;
    uint32_t lastbeg = 0, lastend = 0;
    if (bed->m == 1) {
	lastbeg = bed->a[0]>>32;
	lastend = (uint32_t)bed->a[0];
	length = lastend - lastbeg;
	goto mark;
    }
    reg_sort(bed);
    uint64_t *b;
    b = (uint64_t*) needmem((bed->m+1) * sizeof(uint64_t));

    for (i = 0; i < bed->m; ++i) {
	uint32_t beg, end;
	beg = bed->a[i]>>32; end = (uint32_t)bed->a[i];
	if (lastend < 1) {
	    lastend = end;
	    lastbeg = beg;
	    continue;
	}
	if (lastend >= beg) {
	    if (lastend < end) lastend = end;
	    bed->count[m]++;
	} else {
	    b[m++] = (uint64_t) lastbeg<<32 |lastend;
	    length += lastend - lastbeg;
	    lastbeg = beg; lastend = end;
	}
    }
    if (lastend > 0) {
	length += lastend - lastbeg;
	b[m++] = (uint64_t) lastbeg<<32 | lastend;
    }
    memset(bed->a, 0, bed->m*sizeof(uint64_t));
    memcpy(bed->a, b, m*sizeof(uint64_t));
    bed->m = m;
    freemem(b);
mark:
    bed->flag = 0;
    bed->sorted |= BD_IS_FMTED;
    bed->l_reg = length;
    bed->l_chr = 0;
}

/*
  return :
  0    find this pos in bed region
  -1   in the out range of the regions
  -2   debug
  int  gap with the nearest block
 */
int pos_find(reglist_t *bed, uint32_t pos)
{
#ifdef _DEBUG_MODE
    debug("[%s]", __func__);
#endif
    
    int low = 0;
    int high = bed->m;
    int mid;
    uint32_t beg = (uint32_t)(bed->a[0]>>32);
    uint32_t end = (uint32_t)bed->a[0];
    if (pos < beg && pos > end) return -1;
    if (pos >= beg && pos <= end) return 0;
    while (low < high)
    {
	mid = (high+low)/2;
	beg = (uint32_t)(bed->a[mid]>>32);
	end = (uint32_t)bed->a[mid];
	if (pos >= beg && pos <= end) return 0;
	if (pos < beg) high = mid -1;
	if (pos > end) low = mid + 1;
    }
    if (low == high) {
	beg = (uint32_t)(bed->a[low]>>32);
	end = (uint32_t)bed->a[low];
	uint32_t lastend = (uint32_t)bed->a[low-1];
	uint32_t nextbeg = (uint32_t)(bed->a[low+1]>>32);
	return pos < beg ? beg -pos > pos -lastend ? pos-lastend : beg-pos : pos - end > nextbeg - pos ? nextbeg - pos : pos - end;
    }
    return -2;
}

static reglist_t * reg_clone(reglist_t * bed)
{
#ifdef _DEBUG_MODE
    debug("[%s]", __func__);
#endif
    assert (bed != NULL);
    //if ( bed == NULL ) return NULL;
    reglist_t *bed1 = (reglist_t*)calloc(1, sizeof(reglist_t));
    memcpy(bed1, bed, sizeof(reglist_t));
    if ( bed->n ) {
       	bed1->a = (uint64_t*)calloc(bed->n, sizeof(uint64_t));
	memcpy(bed1->a, bed->a, bed->n *sizeof(uint64_t));
    }
    // index and data
    return bed1;
}

static reglist_t * reg_merge(reglist_t ** beds, int n_beds)
{
#ifdef _DEBUG_MODE
    debug("[%s] n_beds : %d", __func__,n_beds);
#endif
    if (beds == NULL || n_beds == 0) return NULL;
    reglist_t *bed = reg_add(beds, n_beds);
    if (bed==NULL) return NULL;
    regcore_merge(bed);
    return bed;
}

static reglist_t * reg_uniq(reglist_t ** regs, int n_bed ) 
{
#ifdef _DEBUG_MODE
    debug("[%s] n_bed : %d", __func__,n_bed);
#endif

    int i, l;
    if ( n_bed < 2 ) return NULL;
    reglist_t * bed = NULL;
    int id = -1;
    for ( l = 0; l < n_bed; ++l )
    {
	if (regs[l] == NULL) {
	    if(bed) destroy_reg(bed);
	    return NULL;
	}
	if ( id  == -1 )
	    id = regs[l]->id;
	else if (id != regs[l]->id )
	    errabort("bed->id != regs[%d]->id",l);
	
	if ( !(regs[l]->sorted & BD_IS_MERGED) )
	    regcore_merge(regs[l]);
	reglist_t * a[] = { bed, regs[l] };
	reglist_t * newbed = reg_add(a, 2);
	if ( bed ) {
	    destroy_reg(bed);
	    bed=NULL;
	    //clear_reg(bed, destroy_void);
	    //freemem(bed);
	} else {
	    bed = newbed;
	    if ( bed->sorted != BD_IS_SORTED ) regcore_merge(bed);
	    continue;
	}
	bed = newbed;
	reg_sort(bed);    
	uint32_t lastbeg = 0;
	uint32_t lastend = 0;
	uint64_t *b;
	b = (uint64_t*)needmem((bed->m+1) * sizeof(uint64_t));
	int j = 0;
	uint32_t length = 0;
	for (i = 0; i < bed->m; ++i) {
	    int beg, end;
	    beg = bed->a[i]>>32; end = (uint32_t)bed->a[i];

	    /* condition 1:  init beg and end
	     *
	     *         lastbeg        lastend
	     *         |              |
	     *         ===============
	     *                            ---------------------
	     *                           |                    |
	     *                           beg                  end
	     *
	     */
	    if (lastend <= beg || lastend < 1) {
		lastbeg = beg; lastend = end;
		continue;
	    }
    
	    /* condition 2:  init beg and end
	     *
	     *         lastbeg                   lastend
	     *         |                              |
	     *         ================================
	     *              ---------------------
	     *             |                    |
	     *             beg                  end
	     *
	     */
    
	    if (lastend > end) {
		length += end - beg;
		b[j++] = (uint64_t) beg <<32 | (uint32_t)end;
		lastbeg = end;
		continue;
	    } else {
		/* condition 3:  init beg and end
		 *
		 *         lastbeg        lastend
		 *         |              |
		 *         ===============
		 *                 ---------------------
		 *                |                    |
		 *                beg                  end
		 *
		 */
		length += lastend - beg;
		b[j++] = (uint64_t) beg <<32 | (uint32_t)lastend;
		lastbeg = lastend;
		lastend = lastend == end ? 0 : end;
	    }
	}
	bed->n = bed->m+1;
	bed->m = j;
	freemem(bed->a);
	bed->a = b;
	bed->l_reg = length;
	bed->sorted |= BD_IS_FMTED;
    }
    //bed->sorted = BD_IS_SORTED;
    bed->id = id;
    return bed;
}

static reglist_t * reg_diff(reglist_t ** regs, int n_regs)
{
#ifdef _DEBUG_MODE
    debug("[%s] n_regs: %d", __func__,n_regs);
#endif

    int i, j = 0;
    assert ( n_regs > 0 );
    if ( n_regs ==  1 ) {
	reglist_t *bed = reg_clone(regs[0]);
	regcore_merge(bed);
	return bed;
    }
    if (regs[0] == NULL) return NULL;
    reglist_t * mrgs = reg_merge( regs+1, n_regs -1 );
    if (mrgs == NULL) 
	return reg_clone(regs[0]);
    
    reglist_t * tmp[] = { regs[0], mrgs };
    reglist_t * uniq = reg_uniq(tmp, 2);
    //regcore_merge(uniq);
    reglist_t * tmp1[] = { regs[0], uniq };
    reglist_t * reg = reg_add(tmp1, 2);
    destroy_reg(mrgs);
    destroy_reg(uniq);
    reg_sort(reg);
    uint64_t *b;
    b = (uint64_t*)needmem((reg->m+1) * sizeof(uint64_t));
    uint32_t lastbeg = 0; uint32_t lastend = 0;
    uint32_t length = 0;
    for (i = 0; i < reg->m; ++i) {
	uint32_t beg, end;
	beg = reg->a[i]>>32; end = (uint32_t)reg->a[i];

	/* condition 1:  init beg and end
	 *
	 *         beg             end
	 *         |              |
	 *         ===============
	 *
	 */
	if (lastend < 1) {
	    lastbeg = beg; lastend = end;
	    continue;
	}

	/* condition 2: 
	 *
	 *         lastbeg     lastend
	 *         |          |
	 *  region -----------
	 *                 (no overlap)
	 *                           ===============  new region
	 *                          |              |
	 *                          beg            end
	 */
	if (lastend <= beg) {
	    length += lastend - lastbeg;
	    b[j++] = (uint64_t) lastbeg << 32| lastend;
	    lastbeg = beg; lastend = end;
	    continue;
	}
	
	if (lastbeg == beg) {      
	    /* condition 3: 
	     *
	     *         lastbeg     lastend
	     *         |               |
	     *  region ----------------
	     *         |||||||||||||||| (equal)
	     *         ===============  new region
	     *        |              |
	     *        beg            end
	     */
	    if (end == lastend) {
		lastbeg = lastend = 0;
		continue;
	    }

	    /* condition 4: 
	     *
	     *         lastbeg     lastend
	     *         |          |
	     *  region -----------
	     *         | overlap |~~~~~~~
	     *         ==================  new region
	     *        |                 |
	     *        beg            end
	     */
      
	    if (lastend < end) {
		//lastbeg = lastend + 1;
		lastbeg = lastend;
		lastend =  end;
	    } else {
		/* condition 5: 
		 *
		 *         lastbeg                    lastend
		 *         |                          |
		 *  region ---------------------------
		 *         | overlap |||||||||
		 *         ==================  new region
		 *         |                |
		 *         beg            end
		 */
		lastbeg = end;
	    }

	    assert(lastbeg < lastend);
	    continue;
	}
	// lastend > beg come here
	if (lastbeg < beg) {
	    //b[j++] = (uint64_t) lastbeg << 32| (beg -1);
	    length += beg - lastbeg;
	    b[j++] = (uint64_t) lastbeg << 32| beg ; // 0based beg
	} else {
	    errabort("FIXME: not properly sorted!"
	     "Contact developer if you see this message!");
	}
    
	/* condition 6: 
	 *
	 *         lastbeg              lastend
	 *         |                   |
	 *  region --------------------
	 *         ***********||||||||~~~~~~~~~~
	 *                    ==================  new region
	 *                   |                |
	 *                   beg            end
	 */
	if (lastend < end) {
	    //lastbeg = lastend + 1;
	    lastbeg = lastend;
	    lastend = end;
	    continue;
	}

	/* condition 7: 
	 *
	 *         lastbeg                 lastend
	 *         |                            |
	 *  region -----------------------------
	 *         ***********||||||||||||||||||
	 *                    ==================  new region
	 *                   |                |
	 *                   beg            end
	 */
	if (lastend == end) {
	    lastbeg = lastend = 0;
	    continue;
	} else if (lastend > end) {
	    /* condition 8: 
	     *
	     *         lastbeg                                lastend
	     *         |                                            |
	     *  region ---------------------------------------------
	     *         ***********||||||||||||||||||~~~~~~~~~~~~~~~~
	     *                    ==================  new region
	     *                   |                |
	     *                   beg            end
	     */
	    lastbeg = end;
	    continue;
	}
    }

    if ( lastend > 0) {     /*  bugs FIX 2015/03/28 */
	length += lastend - lastbeg;
	b[j++] = (uint64_t) lastbeg << 32| lastend ; // 0based beg
    }
    reg->m = j;
    freemem(reg->a);
    reg->a = b;
    reg->sorted |= BD_IS_FMTED;
    reg->l_reg = length;
    return reg;
}

static reglist_t * reg_comp(reglist_t **regs, int n_regs)
{
#ifdef _DEBUG_MODE
    debug("[%s] n_regs: %d", __func__, n_regs);
#endif
    
    if ( n_regs == 1 ) {
	reglist_t *bed = reg_clone(regs[0]);
	regcore_merge(bed);
	return bed;	
    }
    reglist_t *m = reg_merge( regs, 1);
    if ( m == NULL ) return NULL;
    reglist_t *mrgs = reg_merge( regs+1, n_regs -1 );
    m->data = (void*)mrgs;
    return m;
}

struct cut_pos {
    uint32_t pos;
    struct cut_pos *next;
};
    
static reglist_t *reg_split(reglist_t **reg, int n_regs)
{
#ifdef _DEBUG_MODE
    debug("[%s] n_reg : %d", __func__,n_regs);
#endif

    reglist_t *rb = reg_add(reg, n_regs);
    reg_sort(rb);
    uint64_t *b;
    int n = rb->n > 0 ? rb->n : 2;
    b = (uint64_t*)needmem(n * sizeof(uint64_t));
    int i;
    int m = 0;
    uint32_t lastbeg = 0;
    struct cut_pos *lastends = NULL;
    for (i = 0; i < rb->m; ++i) {
	uint32_t beg = rb->a[i]>>32;
	uint32_t end = (uint32_t)rb->a[i];
	struct cut_pos *curr_end = (struct cut_pos*)malloc(sizeof(struct cut_pos));
	curr_end->pos = end;
	curr_end->next = NULL;
	if (lastends == NULL) {
	    if (lastbeg == 0) {
		lastbeg = beg;
		lastends = curr_end;
		continue;
	    } else {
		errabort("FIXME: lastends is NULL");
	    }
	}
	while (lastbeg < beg) {
	    if ( lastends == NULL ) {
		lastbeg = beg;
		continue;
	    }
	    struct cut_pos * tmp = lastends;
	    uint32_t cut = tmp->pos > beg ? beg : tmp->pos;
	    while (tmp) {
		if (beg < cut) {
		    b[m++] = (uint64_t)lastbeg<<32|beg;
		} else {
		    b[m++] = (uint64_t)lastbeg<<32|cut;
		}
		tmp = tmp->next;
		if ( m == n ) {
		    n = m*2;
		    b = (uint64_t*)realloc(b, n*sizeof(uint64_t));
		}
	    }
	    lastbeg = cut;
	    while (lastends && cut == lastends->pos) {
		struct cut_pos *tmp1 = lastends;
		lastends = lastends->next;
		freemem(tmp1);
	    }
	    if (lastends == NULL) lastbeg = beg;
	}
	if (lastends) {
	    struct cut_pos *tmp = lastends;
	    while (tmp->next && tmp->pos < end) tmp = tmp->next;
	    if (lastends->pos > end) {
		curr_end->next = lastends;
		lastends = curr_end;
	    } else {
		curr_end->next = tmp->next;
		tmp->next = curr_end;
	    }
	} else {
	    lastends = curr_end;
	}
    }

    while (lastbeg && lastends) {
	struct cut_pos *tmp = lastends;
	uint32_t cut = tmp->pos;
	while(tmp)
	{
	    if (lastbeg < cut ) {
		b[m++] = (uint64_t)lastbeg<<32|cut;
	    } else {
		errabort("FIXME: %d\t%d", __func__, __LINE__, lastbeg, cut);
	    }
	    if ( m == n ) {
		n = m*2;
		b = (uint64_t*)realloc(b, n*sizeof(uint64_t));
	    }
	    tmp = tmp->next;
	}
	lastbeg = lastends->pos;
	while (	lastends && lastbeg == lastends->pos)
	{
	    tmp = lastends;
	    lastends = lastends->next;
	    freemem(tmp);
	}
    }
    freemem(rb->a);
    rb->a = b;
    rb->m = m;
    rb->n = n;
    rb->flag ^= BD_IS_SORTED;
    reg_sort(rb);
    return rb;
}

typedef reglist_t * (*handle_func)(reglist_t ** regs, int n_regs);

static bedaux_t * bed_handle(bedaux_t *beds, handle_func func, int *check_error)
{
#ifdef _DEBUG_MODE
    debug("[%s] func: %p", __func__, func);
#endif
    
    bedaux_t * beda = (bedaux_t *)calloc(1, sizeof(bedaux_t));
    beda->alloced = 2;
    beda->is_empty = BD_IS_EMPTY;
    beda->n_files = 1;
    beda->n_seq = beds->n_seq;
    //*beda = { BD_IS_HAVE, 0, NULL, 1, NULL, NULL, beds->n_seq, 0, 0 };
    beda->filenames = (char **)calloc(beda->alloced, sizeof(char*));
    beda->hfiles = (bedfile_t **)calloc(beda->alloced, sizeof(bedfile_t*));
    beda->seq_names = (char**)calloc(beda->n_seq, sizeof(char*));
    beda->hfiles[0] = (bedfile_t*)calloc(1, sizeof(bedfile_t));
    beda->hfiles[0]->reg = kh_init(reg);
    beda->region = beda->length = 0;
    int i,j;
    //int skip_tag;
    for ( i = 0; i < beds->n_seq; ++i )
    {
	reglist_t *regs[beds->n_files];
	beda->seq_names[i] = strdup(beds->seq_names[i]);
	int m = 0, n = 0;
	for ( j = 0; j < beds->n_files; ++j )
	{
	    if ( !beds->hfiles[j] ) {
		regs[n++] = NULL;
		continue;
	    }
	    //if ( beds->hfiles[i]->is_empty == BD_IS_EMPTY ) continue;
	    khiter_t k;
	    regHash_t *rgh = beds->hfiles[j]->reg;
	    k = kh_get(reg, rgh, beds->seq_names[i]);
	    if ( k == kh_end(rgh) ) {
		regs[n++] = NULL;
	    } else {
		regs[n++] = kh_val(rgh, k);
		m++;
	    }
	}
	if ( m == 0 ) continue;
#ifdef _DEBUG_MODE
	debug("[%s] name: %s =======", __func__, beds->seq_names[i]);
#endif
	reglist_t * bed1 = func(regs, n);
	if ( bed1 == NULL ) continue;
	beda->region += bed1->m;
	beda->length += bed1->l_reg;
	if ( *check_error > BD_CHECK_NO && check_reg_by_chromosome_length(bed1) == 1 ) *check_error = BD_CHECK_YES;
	int ret;
	khiter_t k;
	k = kh_put(reg, beda->hfiles[0]->reg, beda->seq_names[i], &ret);
	kh_val(beda->hfiles[0]->reg, k) = bed1;
	beda->is_empty = BD_IS_HAVE;
    }

    return beda;
    
}


static bedaux_t * bed_merge(bedaux_t *beds, int *check_error)
{
#ifdef _DEBUG_MODE
    debug("[%s]", __func__);
#endif
    
    return bed_handle(beds, reg_merge, check_error);
}

static bedaux_t * bed_uniq(bedaux_t *beds)
{
#ifdef _DEBUG_MODE
    debug("[%s]", __func__);
#endif

    int check = BD_CHECK_NO;
    return bed_handle(beds, reg_uniq, &check);
}

static bedaux_t * bed_diff(bedaux_t *beds)
{
#ifdef _DEBUG_MODE
    debug("[%s]", __func__);
#endif

    int check = BD_CHECK_NO;
    return bed_handle(beds, reg_diff, &check);
}

static bedaux_t * bed_count(bedaux_t *beds)
{
#ifdef _DEBUG_MODE
    debug("[%s]", __func__);
#endif
    
    int check = BD_CHECK_NO;
    return bed_handle(beds, reg_split, &check);
}

static void bed_save(const char *fn, bedaux_t * bed) 
{
#ifdef _DEBUG_MODE
    debug("[%s]", __func__);
#endif
    
    FILE *fp;
    fp = open_wfile(fn);
    khiter_t k;
    int i, j;
    bedfile_t *beds = bed->hfiles[0];
    if ( beds == NULL ) {
	fclose(fp);
	return;
    }
    for ( i = 0; i < bed->n_seq; ++i )
    {
	k = kh_get(reg, beds->reg, bed->seq_names[i]);
	if ( k != kh_end(beds->reg)) {
	    reglist_t *reg = kh_val(beds->reg, k);
	    if ( reg == NULL ) continue;
	    for (j = 0; j < reg->m; ++j)
	    {
		uint32_t beg = reg->a[j] >>32;
		uint32_t end = (uint32_t)reg->a[j];
		fprintf(fp, "%s\t%u\t%u\n", bed->seq_names[i], beg, end);
	    }
	}
    }
    fclose(fp);
}


static bedaux_t * bed_comp(bedaux_t * bed)
{
#ifdef _DEBUG_MODE
    debug("[%s]", __func__);
#endif

    int check = BD_CHECK_NO;
    return bed_handle(bed, reg_comp, &check);
}


static bedHandle_t defaultBedHandler =
{
    bed_read,
    bed_merge,
    bed_uniq,
    bed_diff,
    bed_comp,
    bed_count,
    bed_save,
    bed_destroy,
    bed_clear
};

bedHandle_t const *bedHand = &defaultBedHandler;
