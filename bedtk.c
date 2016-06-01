#include "bedutil.h"
#include "commons.h"

const char * version = "0.1.3";

extern bedHandle_t *bedHand;

static bool one_based = FALSE;

static int flank_right = 0;
static int flank_left = 0;
static int trim_tag = BD_IS_FLANK;
static int help = 0;
static char *out = NULL;

int init_argv(int argc, char **argv, bedaux_t *bed)
{
    int n, i;
    while ((n = getopt(argc, argv, "o:hr:l:1")) >= 0)
    {
	switch(n) {	
	    case 'o': out = optarg; break;
	    case 'h': help = 1; break;
	    case 'r': flank_right = atoi(optarg); break;
	    case 'l': flank_left = atoi(optarg); break;
	    case '1': one_based = TRUE; break;
	    default: errabort("invaild option -- %c", n);
	}
    }
    if ( help ) return help;
    n = argc - optind;
    if (n == 0) return 1;
    int ret = 0;
    for (i = 0; i < n; ++i) {
	if ( one_based ) flank_left += 1;
	bedHand->read(argv[optind+i], bed, flank_right, flank_left, &ret, trim_tag);
    }
    if (!one_based && ret) {
	warnings("This file might be not a standard bed format."
		 "Please use parameter \"-1\" if the start in the region file is 1-based!");
    }
    return 0;
}

int mergeHelp(int tag)
{
    if ( tag == BD_IS_TRIM )
	fprintf(stderr, "Usage : program   trim   <in1.bed>\n");
    else
	fprintf(stderr, "Usage : program   merge   <in1.bed> [in2.bed] ...\n");
    fprintf(stderr, "==================================================\n"
	    "  -r     add n bp before each regions\n"
	    "  -l     add n bp after each regions\n"
	    "  -o     output file\n"
	    "  -1     the input bed file is 1-based\n"
	    "  -h     show this message\n"
	);
    return 1;
}

int mergeBed(int argc, char * argv[], int trim_tag) 
{
    bedaux_t bed = INIT_BED_EMPTY; 
    if( init_argv(argc, argv, &bed) ) return mergeHelp(trim_tag);
    int check = BD_CHECK_NO;
    bedaux_t * bed1 = bedHand->merge(&bed, &check);
    if ( out ) bedHand->save(out, bed1);
    assert(bed.region >= bed1->region);
    uint32_t merge_regions = bed.region - bed1->region;
    fprintf(stderr, "Merged %u regions.\n"
	    "Total number of regions : %u.\n"
	    "Length of these regions : %u bp.\n",
	    merge_regions, bed1->region, bed1->length);
    bedHand->destroy(bed1, destroy_void);
    bedHand->clear(&bed, destroy_void);
    return 0;
}

int uniqHelp()
{
    fprintf(stderr,"Usage: uniq <in1.bed> <in2.bed>\n"
	    "   -o  Set output.\n"
	    "   -h  See this information.\n"
	);
    return 1;
}

int uniqBed(int argc, char * argv[])
{
    trim_tag = BD_IS_FLANK;
    bedaux_t bed = INIT_BED_EMPTY;
    if( init_argv(argc, argv, &bed) ) return uniqHelp();
    if ( bed.n_files < 2 )
    {
	bedHand->clear( &bed, destroy_void );
	errabort("At least 2 files for uniq.");    
    }
  
    bedaux_t *bed1 = bedHand->uniq(&bed);
    if ( out ) bedHand->save(out, bed1);
    fprintf(stderr, "Total number of regions : %u.\n"
	    "Length of the uniq regions : %d bp.\n",
	    bed1->region, bed1->length);
    bedHand->destroy( bed1, destroy_void );
    bedHand->clear( &bed, destroy_void );
    return 1;
}

int diffHelp()
{
    fprintf(stderr,"Usage: diff <in1.bed> <in2.bed>\n"
	    "   -o  Set output file.\n"
	    "   -h  See this message.\n"
	);
    return 1;
}

int diffBed(int argc, char * argv[]) 
{
    trim_tag = BD_IS_FLANK;
    bedaux_t bed = INIT_BED_EMPTY;
    if( init_argv(argc, argv, &bed) ) return uniqHelp();
    if ( bed.n_files < 2 )
    {
	bedHand->clear( &bed, destroy_void );
	errabort("At least 2 files for uniq.");    
    }  
    bedaux_t *bed1 = bedHand->diff(&bed);
    if ( out ) bedHand->save(out, bed1);
    fprintf(stderr, "Total number of regions : %u.\n"
	    "There are %d bp included in the first bed files compared with others.\n",
	    bed1->region, bed1->length);
    bedHand->destroy( bed1, destroy_void );
    bedHand->clear( &bed, destroy_void );
    return 0;
}


int trimBed(int argc, char *argv[])
{
    trim_tag = BD_IS_TRIM;
    bedaux_t bed = INIT_BED_EMPTY;
    if( init_argv(argc, argv, &bed) ) return mergeHelp(trim_tag);
    int check = BD_CHECK_NO;
    bedaux_t * bed1 = bedHand->merge(&bed, &check);
    if ( out ) bedHand->save(out, bed1);
    assert(bed.length >= bed1->length);
    uint32_t trimmed_length = bed.length - bed1->length;

    fprintf(stderr, "Trimmed %u bp.\n"
	    "Total number of regions : %u.\n"
	    "Length of these regions : %u bp.\n",
	    trimmed_length, bed1->region, bed1->length);
    bedHand->destroy(bed1, destroy_void);
    bedHand->clear(&bed, destroy_void);
    return 0;
}

int compHelp()
{
    fprintf(stderr,"Usage: comp <in1.bed> <in2.bed>\n");
    return 1;
}

void compare_regions(bedaux_t *bed)
{
    if ( bed == NULL ) return;
    int i, j, l;
    regHash_t *rgh = bed->hfiles[0]->reg;
    writeout("#chr\tstart\tend\tlength\tcovered\tcover_rate\tcover_start\tcover_end\tcount\n");
    for ( l = 0; l < bed->n_seq; ++l )
    {
	khiter_t k;
	k = kh_get(reg, rgh, bed->seq_names[l]);
	if ( k == kh_end(rgh) ) continue;
	reglist_t *reg = kh_val(rgh, k);
	reglist_t *comp = (reglist_t*)reg->data;

	/* chr\tstart\tend\t\length\tcover_length\tcover_rate\tcover_start\t%cover_end */
	if ( comp == NULL ) {
	    for ( i = 0; i < reg->m; ++i )
	    {
		uint32_t start =(uint32_t)(reg->a[i] >> 32);
		uint32_t end = (uint32_t)reg->a[i];
		writeout("%s\t%u\t%u\t%u\t0\t0\tNA\tNA\t0\n", bed->seq_names[l], start, end, end - start);
	    }
	} else {
	    i = j = 0;
	    bool is_uncover = TRUE;
	    while ( i < reg->m && j < comp->m )
	    {
		uint32_t count = comp->count[j];
		uint32_t start =(uint32_t)(reg->a[i] >> 32);
		uint32_t end = (uint32_t)reg->a[i];
		uint32_t cstart =(uint32_t)(comp->a[j] >> 32);
		uint32_t cend = (uint32_t)comp->a[j];
		/*
		  sitution 1 :
		  cstart              cend
		  ====================
		                           ------------------
					   start             end
		 */
		if ( cend <= start ) {
		    j++;
		    continue;
		}
		/*
		  sitution 2 :
		                          cstart              cend
					  ====================
		  ------------------
		  start             end
		*/
		if ( end <= cstart ) {
		    i++;
		    if ( is_uncover )writeout("%s\t%u\t%u\t%u\t0\t0\tNA\tNA\t0\n", bed->seq_names[l], start, end, end - start);
		    is_uncover = TRUE;
		    continue;
		}
		/*
		  sitution 3 :
		  cstart              cend
		  ====================
		                    ------------------
	      			   start             end

		  sitution 5 :
		       cstart        cend
		           ==========
	               ------------------
	   	   start             end

		 */
		if ( cend > start && cend < end) {
		    j++;
		    is_uncover = FALSE;
		    uint32_t length = end - start;
		    uint32_t cov = cstart < start ? cend - start : cend - cstart;
		    writeout("%s\t%u\t%u\t%u\t%u\t%.6f\t%u\t%u\t%u\n", bed->seq_names[l], start, end, length, cov, (float)cov/length, cstart, cend, count);
		    continue;
		}
		/*
		  sitution 4 :
		  cstart                                    cend
		  ===========================================
		                    ------------------
	      			   start             end

		  sitution 6 :
		                          cstart              cend
			      	      ===================
	               ------------------
	   	   start             end
		 */

		if ( cend >= end ) {
		    i++;
		    is_uncover = TRUE;
		    uint32_t length = end - start;
		    uint32_t cov = cstart > start ? end - cstart : end - start;
		    writeout("%s\t%u\t%u\t%u\t%u\t%.6f\t%u\t%u\t%u\n", bed->seq_names[l], start, end, length, cov, (float)cov/length, cstart, cend, count);
		    continue;
		}
		debug("Never comes here! %s\t%u\t%u\t%u\t%u\n", bed->seq_names[l], start, end, cstart, cend);
	    }
	    while ( i < reg->m )
	    {
		uint32_t start =(uint32_t)(reg->a[i] >> 32);
		uint32_t end = (uint32_t)reg->a[i];
		if ( is_uncover ) writeout("%s\t%u\t%u\t%u\t0\t0\tNA\tNA\t0\n", bed->seq_names[l], start, end, end - start);
		is_uncover = TRUE;
		i++;
	    }
	}
    }
}
int compBed(int argc, char *argv[])
{
    trim_tag = BD_IS_FLANK;
    bedaux_t bed = INIT_BED_EMPTY;
    if( init_argv(argc, argv, &bed) ) return compHelp();
    bedaux_t * bed1 = bedHand->comp(&bed);
    compare_regions(bed1);
    bedHand->destroy(bed1, destroy_reg);
    bedHand->clear(&bed, destroy_void);
    return 0;
}

int summaryHelp()
{
    fprintf(stderr,"Usage: sum <in.bed>\n");
    return 1;
}

int summaryBed(int argc, char *argv[])
{
    int i;
    bedaux_t bed = INIT_BED_EMPTY;
    int check = BD_CHECK_NO;
    if (init_argv(argc, argv, &bed)) return summaryHelp();

    bedaux_t *bed1 = bedHand->merge(&bed, &check);
    if (bed.n_files > 1) {
	warnings("only summary the first file: %d\n", bed.n_files);
    }
    bedfile_t *b = bed1->hfiles[0];
    regHash_t *rh = b->reg;
    writeout("#chr\tcount\tlength\tpercent\n");
    for (i=0; i<bed1->n_seq; ++i) {
	khiter_t k;
	k = kh_get(reg, rh, bed1->seq_names[i]);
	if (k == kh_end(rh))
	    continue;
	reglist_t *r = kh_val(rh, k);
	if (r) {
	    writeout("%s\t%u\t%u\t%.4f\n", bed1->seq_names[i], r->m, r->l_reg, (float)r->l_reg/bed1->length);
	}
    }

    bedHand->destroy(bed1, destroy_reg);
    bedHand->clear(&bed, destroy_void);
    return 0;
}

int lengthHelp()
{
    fprintf(stderr,"Usage: sum <in.bed>\n");
    return 1;
}

int lengthBed(int argc, char *argv[])
{
    int i,j;
    bedaux_t bed = INIT_BED_EMPTY;
    int check = BD_CHECK_NO;
    if (init_argv(argc, argv, &bed)) return lengthHelp();

    bedaux_t *bed1 = bedHand->merge(&bed, &check);

    bedfile_t *b = bed1->hfiles[0];
    regHash_t *rh = b->reg;
    writeout("#chr\tstart\tend\tlength\tpchr\tgap_length\n");
    for (i=0; i<bed1->n_seq; ++i) {
	khiter_t k;
	k = kh_get(reg, rh, bed1->seq_names[i]);
	if (k == kh_end(rh))
	    continue;
	reglist_t *r = kh_val(rh, k);
	if (r) {
	    uint32_t lastend = 0;
	    for (j=0; j<r->m; ++j) {
		uint32_t begin = (uint32_t)(r->a[j]>>32);
		uint32_t end = (uint32_t)(r->a[j]);
		uint32_t gap = lastend == 0 ? 0 : begin - lastend;
		lastend = end;
		uint32_t length = end - begin;
		writeout("%s\t%u\t%u\t%u\t%.4f\t%u\n", bed1->seq_names[i], begin, end, length, (float)length/r->l_reg, gap);
	    }
	}
    }

    bedHand->destroy(bed1, destroy_reg);
    bedHand->clear(&bed, destroy_void);
    return 0;
    
}
int countHelp()
{
    fprintf(stderr,"Usage: count <in1.bed> <in2.bed>\n");
    return 1;
}

void count_regions(bedaux_t *bed)
{
    if ( bed == NULL ) return;
    int i, l;
    regHash_t *rgh = bed->hfiles[0]->reg;
    writeout("#chr\tstart\tend\tcount\n");
    for ( l = 0; l < bed->n_seq; ++l )
    {
	khiter_t k;
	k = kh_get(reg, rgh, bed->seq_names[l]);
	if ( k == kh_end(rgh) ) continue;
	reglist_t *reg = kh_val(rgh, k);
	uint32_t lastbeg, lastend = 0;
	uint32_t count = 0;
	for ( i= 0; i < reg->m; ++i)
	{
	    uint32_t beg = reg->a[i]>>32;
	    uint32_t end = (uint32_t)reg->a[i];

	    if (lastend == 0) {
		lastbeg = beg;
		lastend = end;
		count++;
		continue;
	    }
	    if (lastbeg == beg && lastend == end) {
		count++;
	    } else if (beg >= lastend) {
		writeout("%s\t%u\t%u\t%u\n", bed->seq_names[l], lastbeg, lastend, count);
		lastbeg = beg;
		lastend = end;
		count = 1;
	    } else {
		errabort("FIXME: lastbeg : %u, lastend : %u, beg : %u, end : %u", lastbeg, lastend, beg, end);
	    }
	}
	if (lastend > 0) {
	    writeout("%s\t%u\t%u\t%u\n", bed->seq_names[l], lastbeg, lastend, count);
	}

    }
    
}
int countBed(int argc, char *argv[])
{
    trim_tag = BD_IS_FLANK;
    bedaux_t bed = INIT_BED_EMPTY;
    if( init_argv(argc, argv, &bed) ) return countHelp();
    bedaux_t * bed1 = bedHand->count(&bed);
    count_regions(bed1);
    bedHand->destroy(bed1, destroy_reg);
    bedHand->clear(&bed, destroy_void);
    return 0;    
}
static int usage(void)
{
    fprintf(stderr,
	    "bedutils -- a toolkit to handle 3 columns bed files, the basic functions of this toolkit are:\n"
	    "=============================================================================================\n"
	    "merge    : merge all the input bed files. use -l and -r to flank the regions.\n"
	    "diff     : find out the different region between first bed file and others.\n"
	    "uniq     : find out the uniq regions of all the input bed files\n"
	    "trim     : trim the input bed files with -r and -l.\n"
	    "comp     : compare the first bed file with others\n"
	    "sum      : summarize the chromosome length in the bed file(s).\n"
	    "length   : summarize the length of each region in the bed file(s).\n"
	    "count    : count the regions \n"
	    "=============================================================================================\n"
	    "Usage:\n"
	    "  merge -o <FILE> <in1.bed [in2.bed ...]>\n"
	    "  diff  -o <FILE> <in1.bed> <in2.bed> [in3.bed ...]\n"
	    "  uniq  -o <FILE> <in1.bed> <in2.bed> [in3.bed ...]\n"
	    "  trim  -l xx -r xx <in1.bed> [in3.bed ...]\n"
	    "  comp  <in.bed> <in2.bed> \n"
	    "  sum   <in.bed>\n"
	    "  length <in.bed>\n"
	    "  count  <in.bed> <in2.bed> \n"
	    "=============================================================================================\n"
	    "Version : %s\n"
	    "Author: Shi Quan (shiquan@genomics.cn)\n"
	    "Website:\n"
	    "https://github.com/shiquan/bedutils\n"
	    , version);
    return 1;
}

int main(int argc, char *argv[]) 
{
    if (argc < 2) return usage();
    else if (STREQ(argv[1], "merge")) return mergeBed(argc-1, argv+1, trim_tag);
    else if (STREQ(argv[1], "diff")) return diffBed(argc-1, argv+1);
    else if (STREQ(argv[1], "uniq")) return uniqBed(argc-1, argv+1);
    else if (STREQ(argv[1], "trim")) return trimBed(argc-1, argv+1);
    else if (STREQ(argv[1], "comp")) return compBed(argc-1, argv+1);
    else if (STREQ(argv[1], "count")) return countBed(argc-1, argv+1);
    else if (STREQ(argv[1], "sum")) return summaryBed(argc-1, argv+1);
    else if (STREQ(argv[1], "length")) return lengthBed(argc-1, argv+1);
    else return usage();
    return 1;
}
