//
//  bedutil.h
//
// Copyright (c) 2014, 2015 Beijing Genomics Institution (BGI)
//
//  Author:      Shi Quan (shiquan@genomics.cn)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
// 
//
// This program is designed to handle bed file, inspired by liheng's bedidx.c
//
#ifndef BEDUTILS_H
#define BEDUTILS_H
#include <stdint.h>
#include "khash.h"

#define BD_CHECK_NO  0 // check the region of each chromosome
#define BD_CHECK_YES 1
#define BD_CHECK_ERR 2

#define BD_MAIN_REG 0 // the main struct of all the bedfile_t[] , usually this first one bedfile_t[0]

#define BD_IS_HAVE  0 // check this bedfile_t is full or empty
#define BD_IS_EMPTY 1

#define BD_IS_UNSORT 0 // check the bed struct is sorted or not when try to sort it
#define BD_IS_SORTED 1
#define BD_IS_MERGED 2
#define BD_IS_FMTED (BD_IS_SORTED|BD_IS_MERGED)

#define BD_IS_FLANK 0
#define BD_IS_TRIM 1

#define maxlen(b) ((b)->m > 0 ? (uint32_t)(b)->a[(b)->m - 1] : 0)

#define BED_MAX_BIN 37450

typedef struct {
    int m; // used-length
    int n; // alloced-length
    int flag;
    int sorted; // sorted flag, check this flag before sort
    int id; // name id
    uint32_t l_reg; // total length of these regions
    uint32_t l_chr; // length of the chromosome, used for check the positions of each region, 0 is not set
    uint64_t *a;
    uint32_t *count; // how many regions merged in this region, for stat only
    int *idx;
    void *data; // self defined data
} bedreglist_t;

typedef bedreglist_t * bedreglist_p;

KHASH_MAP_INIT_STR(reg, bedreglist_p)

typedef kh_reg_t regHash_t;

typedef struct {
    int is_empty;
    uint32_t region;
    uint32_t length;
    regHash_t *reg;
} bedfile_t;

typedef struct {
    int is_empty;
    int alloced;
    int n_files; // file numbers, length of **reg
    char **filenames;
    bedfile_t **hfiles; // structs of bed files
    char **seq_names; // chromosome names
    int n_seq; // number of names
    uint32_t region; // region count
    uint32_t length; // total length of regions
} bedaux_t;

typedef void (*bedvoid_destroy)(void *data);
void destroy_void(void *data);
void destroy_reg(void *data);
//typedef void (*handle_func)(regHash_t *bed);

struct _bedHandle {
    /* handle whole struct */ 
    void (*read)(const char *fn, bedaux_t *bed, int a, int b, int *ret, int trim_tag);
    bedaux_t *(*merge)(bedaux_t *bed, int *is_error);
    bedaux_t *(*uniq)(bedaux_t *bed);
    bedaux_t *(*diff)(bedaux_t *bed);
    bedaux_t *(*comp)(bedaux_t *bed);
    bedaux_t *(*count)(bedaux_t *bed);
    void (*save)(const char *fn, bedaux_t *bed);
    /* destroy */
    void (*destroy)(bedaux_t *bed, bedvoid_destroy func);
    void (*clear)(bedaux_t *bed, bedvoid_destroy func);
};

typedef struct _bedHandle bedHandle_t;

#endif
