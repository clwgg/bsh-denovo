#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>

#include "htslib/htslib/sam.h"
#include "htslib/htslib/khash.h"

typedef struct hdata {
  int sum;
  int pos[4];
} hdata;

KHASH_MAP_INIT_STR(rg, hdata)


typedef struct aux_t {
  samFile *in;
  bam_hdr_t *h;
  hts_idx_t *idx;
  hts_itr_t *itr;
  int min_mapQ;
} aux_t;

typedef struct outfh {
  FILE *map;
  FILE *ped;
  FILE *debug;
} outfh;

typedef struct pos_filt {
  double mis;
  double maj;
  double frq;
} pos_filt;

typedef struct opt_flag {
  int d;
  int s;
  int f;
} opt_flag;

typedef struct lnode {
  char *bases;
  struct lnode *next;
} lnode;

typedef struct llist {
  int npos;
  int n;
  lnode *start;
  lnode *end;
  const char **names;
} llist;

void init_random(void)
{
  char *env = getenv("RANDOM_SEED");
  if (env == NULL || strcmp(env, "") == 0) {
    srand(time(NULL));
    srand48(time(NULL));
  } else {
    srand(atoi(env));
    srand48(atoi(env));
  }
}

int read_bam(void *data, bam1_t *b)
{

  aux_t *aux = (aux_t*)data;
  int ret;

  while (1) {
    ret = aux->itr ? sam_itr_next(aux->in, aux->itr, b) : sam_read1(aux->in, aux->h, b);
    if (ret < 0) break;
    if ( (int)b->core.qual < aux->min_mapQ ) continue;
    if ( b->core.flag & (BAM_FUNMAP | BAM_FDUP) ) continue;
    break;
  }
  return ret;

}

void free_aux(aux_t data)
{
  if(data.h) bam_hdr_destroy(data.h);
  if(data.idx) hts_idx_destroy(data.idx);
  if(data.itr) hts_itr_destroy(data.itr);
  if(data.in) sam_close(data.in);
}

int read_pos(int n, const bam_pileup1_t *plp, khash_t(rg) *rg_h)
{
  int i;
  khint_t k;

  for (i = 0; i < n; ++i) {

    const bam_pileup1_t *p = plp + i;
    uint8_t *seq = bam_get_seq(p->b);
    uint8_t nuc = bam_seqi(seq, p->qpos);
    uint8_t *rg = bam_aux_get(p->b, "RG");

    k = kh_get(rg, rg_h, (char *)rg+1);

    if (kh_exist(rg_h, k)) {
      switch (nuc) {
      case 1:  kh_value(rg_h, k).pos[0]++ ; kh_value(rg_h, k).sum++ ; break;
      case 2:  kh_value(rg_h, k).pos[1]++ ; kh_value(rg_h, k).sum++ ; break;
      case 4:  kh_value(rg_h, k).pos[2]++ ; kh_value(rg_h, k).sum++ ; break;
      case 8:  kh_value(rg_h, k).pos[3]++ ; kh_value(rg_h, k).sum++ ; break;
      }
    }
  }

  //  printf("%s\t%d\t%d\t%d\t%d\t%d\n", kh_key(rg_h, k), kh_value(rg_h, k).pos[0], kh_value(rg_h, k).pos[1], kh_value(rg_h, k).pos[2], kh_value(rg_h, k).pos[3], kh_value(rg_h, k).sum);
  return 0;
}

int reset_pos(khash_t(rg) *rg_h)
{
  khint_t k;
  for (k = kh_begin(rg_h); k != kh_end(rg_h); ++k) {
    if (kh_exist(rg_h, k)) {
      memset(kh_value(rg_h, k).pos, 0, sizeof kh_value(rg_h, k).pos[0] * 4);
      kh_value(rg_h, k).sum = 0;
    }
  }
  return 0;
}

int init_llist(llist *list, khash_t(rg) *rg_h)
{
  list->n = kh_size(rg_h);

  list->names = malloc(sizeof *list->names * list->n);

  int i = 0;
  khint_t k;
  for (k = kh_begin(rg_h); k != kh_end(rg_h); ++k) {
    if (kh_exist(rg_h, k)) {
      list->names[i] = kh_key(rg_h, k);
      ++i;
    }
  }

  list->npos = 0;
  list->start = NULL;
  list->end = NULL;
  return 0;
}

void free_llist(llist *list)
{
  lnode *node = list->start;
  while (node) {
    free(node->bases);
    lnode *tmp = node;
    node = node->next;
    free(tmp);
  }
  free(list->names);
}

int add_lnode(llist *list, char *bases)
{
  lnode *node = NULL;
  node = malloc(sizeof *node);
  node->bases = bases;
  node->next = NULL;
  if (list->start == NULL) {
    list->start = node;
    list->end = node;
  } else {
    list->end->next = node;
    list->end = node;
  }
  list->npos++;
  return 0;
}

int max_id4(int *pos)
{
  int i;
  int j = -1;
  int max = 0;
  for (i = 0; i < 4; ++i) {
    if (pos[i] > max) {
      max = pos[i];
      j = i;
    }
    else if (pos[i] == max && max != 0) {
      if (rand() % 2) {
        max = pos[i];
        j = i;
      }
    }
  }
  return j;
}

int num_zero4(int *pos)
{
  int i;
  int j = 0;
  for (i = 0; i < 4; ++i) {
    if (pos[i] == 0)
      j++;
  }
  return j;
}

int sum4(int *pos)
{
  int i;
  int sum = 0;
  for (i = 0; i < 4; ++i) {
    sum += pos[i];
  }
  return sum;
}

void netw_sort4(int *pos)
{
#define min(x, y) (x<y?x:y)
#define max(x, y) (x<y?y:x)
#define SWAP(x,y) { const int a = min(pos[x], pos[y]);  \
                    const int b = max(pos[x], pos[y]);  \
                    pos[x] = b; pos[y] = a; }
  SWAP(1, 3);
  SWAP(0, 2);
  SWAP(0, 1);
  SWAP(2, 3);
  SWAP(1, 2);
#undef SWAP
#undef min
#undef max
}

int apply_filt(khash_t(rg) *rg_h, pos_filt *filt)
{
  int j;
  int n = 0;
  int major[] = {0,0,0,0};
  int rg_n = kh_size(rg_h);

  khint_t k;
  for (k = kh_begin(rg_h); k != kh_end(rg_h); ++k) {
    if (kh_exist(rg_h, k)) {
      int sum = kh_value(rg_h, k).sum;
      if (sum == 0) continue;

      if (filt->maj == 0) {
        j = max_id4(kh_value(rg_h, k).pos);
        if (j >= 0) {
          major[j]++;
          n++;
        }
      } else {
        int z = 0;
        for (j = 0; j < 4; ++j) {
          if (kh_value(rg_h, k).pos[j] == 0) continue;
          if (kh_value(rg_h, k).pos[j]/(double)sum >= filt->maj) {
            major[j]++;
            z++;
          }
        }
        if (z) n++;
      }
    }
  }

  if (n <= filt->mis * (double)rg_n)
    return 0;
  if (num_zero4(major) > 2)
    return 0;
  netw_sort4(major);
  if (major[1]/(double)n < filt->frq)
    return 0;
  return 1;
}

int rnd1_4(int *pos, int sum)
{
  int i;

  if (sum < 1)
    return -1;

  double prob = 0.0;
  double rnd = drand48();

  for (i = 0; i < 4; ++i) {
    if (pos[i] == 0)
      continue;
    prob += pos[i]/(double)sum;
    if (prob - rnd >= 0)
      break;
  }

  return i;
}

char ntob(int i)
{
  switch (i) {
  case 0:  return 'A';
  case 1:  return 'C';
  case 2:  return 'G';
  case 3:  return 'T';
  }
  return '0';
}

int bton(char c)
{
  switch (c) {
  case 'A':  return 0;
  case 'C':  return 1;
  case 'G':  return 2;
  case 'T':  return 3;
  }
  return -1;
}

int rnd2_4(int *pos, int sum)
{
  int c, s, i, j, n = 0;

  if (sum < 1)
    return -1;

  int *tmp;
  tmp = malloc(sizeof *tmp * sum);

  for (i = 0; i < 4; ++i) {
    for (j = 0; j < pos[i]; ++j) {
      tmp[n] = i;
      ++n;
    }
  }

  s = rand() % sum;
  c = tmp[s];
  free(tmp);
  return c;
}

char sample_pos(int *pos, int flag)
{
  int s, i;
  int sum = sum4(pos);

  if (flag && sum < 2)
    return '0';
  if (sum == 0)
    return '0';

  if (flag) {
    int count[4] = {0,0,0,0};
    int tmp[4] = {0,0,0,0};
    int z;
    memcpy(tmp, pos, 4 * sizeof *tmp);
    for (i = 0; i < 3 && i < sum; ++i) {
      z = rnd1_4(tmp, sum-i);
      tmp[z]--;
      count[z]++;
    }
    s = max_id4(count);
    if (count[s] < 2) {
      return '0';
    }
  } else {
    s = rnd1_4(pos, sum);
  }

  return ntob(s);
}

int produce_ped(llist *list, outfh out)
{
  int i;
  for (i = 0; i < list->n; ++i) {
    fprintf(out.ped, "%s %s 0 0 0 -9", list->names[i], list->names[i]);
    lnode *node = list->start;
    while (node) {
      fprintf(out.ped, " %c %c", node->bases[i], node->bases[i]);
      node = node->next;
    }
    fprintf(out.ped, "\n");
  }

  return 0;
}

int output_gen(khash_t(rg) *rg_h, int tid, int p, aux_t data, llist *list, outfh out, opt_flag *opt)
{
  fprintf(out.map, "%s\t%s_%d\t0\t%d\n", data.h->target_name[tid], data.h->target_name[tid], p, p);

  char *bases = NULL;
  int rg_n = kh_size(rg_h);

  bases = malloc(sizeof *bases * rg_n);

  int i = 0;
  khint_t k;
  for (k = kh_begin(rg_h); k != kh_end(rg_h); ++k) {
    if (kh_exist(rg_h, k)) {
      bases[i] = sample_pos(kh_value(rg_h, k).pos, opt->s);

      if (out.debug) {
        fprintf(out.debug, "%s\t%d\t%s\t%d\t%d\t%d\t%d\t%c\n", data.h->target_name[tid], p, kh_key(rg_h, k),
                kh_value(rg_h, k).pos[0], kh_value(rg_h, k).pos[1], kh_value(rg_h, k).pos[2], kh_value(rg_h, k).pos[3], bases[i]);
      }
      ++i;
    }
  }

  add_lnode(list, bases);

  return 0;
}

int rg_header(aux_t data, khash_t(rg) *rg_h)
{

  char *buffer = data.h->text;
  char *line;
  int n = 0;

  int ret;
  khint_t k;

  buffer = strstr(buffer, "@RG");
  line = strtok(buffer, "\n");

  while (line) {
    if (line[0] == '@' && line[1] == 'R' && line[2] == 'G') {
      line = strstr(line, "ID:");
      int i = 0;
      char *id = NULL;
      while (line[i] != '\t')
        i++;
      line[i] = '\0';
      id = calloc(i, sizeof *id);
      memcpy(id, line + 3, i - 3);
      k = kh_put(rg, rg_h, id, &ret);
      n++;
    }
    line = strtok(0, "\n");
  }

  for (k = kh_begin(rg_h); k != kh_end(rg_h); ++k) {
    if (kh_exist(rg_h, k)) {
      //      printf("Key: %s\n", kh_key(rg_h, k));
      memset(kh_value(rg_h, k).pos, 0, 4 * sizeof *kh_value(rg_h, k).pos);
      kh_value(rg_h, k).sum = 0;
    }
  }

  return 0;
}

void free_rgh(khash_t(rg) *rg_h)
{
  khint_t k;
  for (k = kh_begin(rg_h); k != kh_end(rg_h); ++k) {
    if (kh_exist(rg_h, k)) {
      free((char*)kh_key(rg_h, k));
    }
  }
}

void close_outf(outfh out)
{
  if (out.map) fclose(out.map);
  if (out.ped) fclose(out.ped);
  if (out.debug) fclose(out.debug);
}

int usage(int r, char **argv)
{
  printf("\nUsage: %s [options] multi_RG.bam\n\n", argv[0]);
  printf("Options:\n");
  printf("\t-d\tProduce stats file for debugging and analytics.\n");
  printf("\t-o\tOutput file base (default: out)\n\n");
  printf("    Position discovery (options here don't apply to the sampling of genotypes):\n");
  printf("\t-m\tData completeness required (per position, across samples, default: 0.5).\n");
  printf("\t-a\tMinimum in-sample frequency a base must have to be considered (within sample, default: highest - random if tied).\n");
  printf("\t-f\tMinimum minor allele frequency across samples (default: 1/n).\n\n");
  printf("    Base sampling:\n");
  printf("\t-c\tUse Consensify method for base sampling.\n");
  printf("\n");
  return r;
}

int main(int argc, char **argv)
{

  int elem;
  opt_flag opt = {0,0,0};
  char *base = "out";
  outfh out = {NULL, NULL, NULL};
  pos_filt filt = {0.5, 0, 0};

  while (( elem = getopt(argc, argv, "o:dcm:a:f:") ) >= 0) {
    switch(elem) {
    case 'o': base = optarg; break;
    case 'd': opt.d = 1; break;
    case 'c': opt.s = 1; break;
    case 'm': filt.mis = atof(optarg); break;
    case 'a': filt.maj = atof(optarg); break;
    case 'f': filt.frq = atof(optarg); break;
    }
  }

  if (argc - optind != 1)
    return usage(1, argv);

  aux_t data = {NULL, NULL, NULL, NULL, 0};

  data.in = sam_open(argv[optind], "r");
  data.h = sam_hdr_read(data.in);
  //  data.idx = sam_index_load(data.in, argv[optind]);
  //  data.min_mapQ = 1;

  if (!data.in || !data.h) {
    free_aux(data);
    printf("\nCannot open in.\n");
    return usage(2, argv);
  }

  char *fn_buf = malloc(sizeof *base * (strlen(base) + 10));
  sprintf(fn_buf, "%s.map", base);
  out.map = fopen(fn_buf, "w");
  sprintf(fn_buf, "%s.ped", base);
  out.ped = fopen(fn_buf, "w");
  if (opt.d) {
    sprintf(fn_buf, "%s.stat", base);
    out.debug = fopen(fn_buf, "w");
    fprintf(out.debug, "SEQ\tPOS\tRG\tA\tC\tG\tT\tSM\n");
  }
  free(fn_buf);

  if (!out.map || !out.ped) {
    free_aux(data);
    close_outf(out);
    printf("\nCannot open out.\n");
    return usage(3, argv);
  }

  init_random();

  khash_t(rg) *rg_h = kh_init(rg);
  rg_header(data, rg_h);

  llist list;
  init_llist(&list, rg_h);

  bam_plp_t iter;
  const bam_pileup1_t *plp;
  int tid, p, n;

  iter = bam_plp_init(read_bam, (void*)&data);
  while ((plp = bam_plp_auto(iter, &tid, &p, &n)) != 0) {
    read_pos(n, plp, rg_h);
    if (apply_filt(rg_h, &filt))
      output_gen(rg_h, tid, p, data, &list, out, &opt);
    reset_pos(rg_h);
  }
  produce_ped(&list, out);
  bam_plp_destroy(iter);

  free_aux(data);
  free_rgh(rg_h); kh_destroy(rg, rg_h);
  free_llist(&list);
  close_outf(out);
  return 0;
}
