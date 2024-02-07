#include "utils.h"
#include "htslib/vcf.h"
#include "htslib/vcfutils.h"
#include "htslib/faidx.h"
#include "gtf.h"
#include "number.h"

int usage()
{
    fprintf(stderr, "TRIBE -gtf gene.gtf in.bcf\n");
    fprintf(stderr, "    -gtf        Annotation file.\n");
    fprintf(stderr, "    -fasta      Genome reference.\n");
    fprintf(stderr, "    -bed        Output A>G editing position.\n");
    fprintf(stderr, "    -out-bcf    Output editing BCF.\n");
    fprintf(stderr, "    -out-fa     Output flanking regions around editing position. For motif enrichment analysis.\n");
    fprintf(stderr, "    -flank      Flank size around editing position. [100]\n");
    return 1;
}
struct args {
    const char *input_fname;
    const char *gtf_fname;
    const char *bed_fname;
    const char *fasta_fname;
    const char *fasta_fname_out;
    const char *bcf_fname;
    
    struct gtf_spec *G;
    FILE *bed_fp_out;
    FILE *fasta_fp_out;
    faidx_t *fai;

    int flank;
} args = {
    .input_fname     = NULL,
    .gtf_fname       = NULL,
    .bed_fname       = NULL,
    .bcf_fname       = NULL,
    .G               = NULL,
    .bed_fp_out      = NULL,
    .fasta_fp_out    = NULL,
    .fasta_fname     = NULL,
    .fasta_fname_out = NULL,
    .fai             = NULL,
    .flank           = 100,
};
int parse_args(int argc, char **argv)
{
    if (argc == 1) return 1;
    const char *flk = NULL;
    
    int i;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;

        // common options
        if (strcmp(a, "-gtf") == 0) var = &args.gtf_fname;
        else if (strcmp(a, "-bed") == 0) var = &args.bed_fname;
        else if (strcmp(a, "-fasta") == 0) var = &args.fasta_fname;
        else if (strcmp(a, "-out-fa") == 0) var = &args.fasta_fname_out;
        else if (strcmp(a, "-out-bcf") == 0) var = &args.bcf_fname;
        else if (strcmp(a, "-flank") == 0) var = &flk;
        else if (strcmp(a, "-h") == 0) return 1;
        
        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }
        if (a[0] == '-' && a[1] != '\0') error("Unknown argument, %s",a);
        if (args.input_fname == NULL) {
            args.input_fname = a;
            continue;
        }
        error("Unknown argument: %s", a);
    }

    if (args.gtf_fname == NULL) error("-gtf is required.");

    args.G = gtf_read_lite(args.gtf_fname);
    if (args.G == NULL) error("Failed to load gtf.");
    if (args.bed_fname != NULL) {
        args.bed_fp_out = fopen(args.bed_fname, "w");
        if (args.bed_fp_out == NULL) error("%s : %s.", args.bed_fname, strerror(errno));
    }

    if (flk != NULL) {
        args.flank = str2int((char*)flk);
        if (args.flank < 0) args.flank *= -1;
    }

    if (args.fasta_fname_out) {
        if (args.fasta_fname == NULL) error("-fasta must be set.");
        args.fai = fai_load(args.fasta_fname);
        if (args.fai == NULL) error("Index %s first!", args.fasta_fname);
        args.fasta_fp_out = fopen(args.fasta_fname_out, "w");
        if (args.fasta_fp_out == NULL) error("%s : %s.", args.fasta_fname_out, strerror(errno));
    }
    return 0;
}
void memory_release()
{
    if (args.bed_fp_out) fclose(args.bed_fp_out);
    if (args.fasta_fp_out) fclose(args.fasta_fp_out);
    if (args.fai) fai_destroy(args.fai);

    gtf_destroy(args.G);
}

static int geno_stat[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
static int geno_stat_ori[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};

static char *reverse_seq(char *s, int l)
{
    char *r = malloc(sizeof(char)*(l+1));
    int i;
    for (i = 0; i < l; ++i ) {
        switch(s[i]) {
            case 'A':
            case 'a':
                r[l-i-1]='T'; break;
            case 'C':
            case 'c':
                r[l-i-1]='G'; break;
            case 'G':
            case 'g':
                r[l-i-1]='C'; break;
            case 'T':
            case 't':
                r[l-i-1]='A'; break;
            case 'N':
            case 'n':
                r[l-i-1]='N'; break;
            default:
                error("Unknown bases, %s",s);
        }
    }
    r[l] = '\0';
    return r;
}

int main(int argc, char **argv)
{
    if (parse_args(argc, argv)) return usage();
    htsFile *fp;
    fp = hts_open(args.input_fname, "r");
    htsFormat type = *hts_get_format(fp);
    if (type.format != vcf && type.format != bcf)
        error("Unsupport input format, only accept BCF/VCF format. %s", args.input_fname);

    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (hdr == NULL) error("Failed to parse header of input.");
    bcf1_t *bcf = bcf_init();
    int *tmpi = NULL;
    int mtmpi = 0;

    htsFile *out = NULL;
    if (args.bcf_fname != NULL) {
        out = hts_open(args.bcf_fname, "wb");
        assert(out);
    }
    
    
    if (out) bcf_hdr_write(out, hdr);
    
    for (;;) {
        if (bcf_read(fp, hdr, bcf)) break;
        if (bcf->rlen != 1) continue; // only consider snp
        if ( bcf_get_variant_types(bcf) == VCF_REF ) continue; // skip ref
        if (bcf->n_allele != 2) continue; // skip mnps

        struct region_itr *itr = gtf_query(args.G, (char*)hdr->id[BCF_DT_CTG][bcf->rid].key, bcf->pos, bcf->pos+1);
        if (itr == NULL || itr->n == 0) continue;
        
        struct gtf *gtf = NULL;
        int i;
        for (i = 0; i < itr->n; ++i) {
            struct gtf *gtf0 = itr->rets[i];
            if (bcf->pos+1 < gtf0->start || bcf->pos+1 > gtf0->end) continue;
            if (gtf != NULL) { 
                gtf = NULL; 
                break; // only keep the first gene or transcript
            }
            gtf = itr->rets[i];
        }

        if (gtf == NULL) continue;
        
        int ref = bcf_acgt2int(*bcf->d.allele[0]);
        int alt = bcf_acgt2int(*bcf->d.allele[1]);
        geno_stat_ori[ref][alt]++;
        
        if (gtf->strand == 1) {
            ref = 3-ref;
            alt = 3-alt;
        }
        geno_stat[ref][alt]++;

        // get the DP value
        int dp = 0;
        bcf_get_info_int32(hdr, bcf, "DP",&tmpi, &mtmpi);
        if (ref == 0 && alt == 2) {
            if (args.bed_fp_out)
                fprintf(args.bed_fp_out, "%s\t%d\t%d\t%s\t%d\t%c\n", (char*)hdr->id[BCF_DT_CTG][bcf->rid].key, (int)bcf->pos, (int)bcf->pos+1,
                        dict_name(args.G->gene_name,gtf->gene_name),mtmpi==1 ? tmpi[0] : 0, gtf->strand == 0 ? '+' : '-');
            if (out) bcf_write(out, hdr, bcf);
            if (args.fasta_fp_out) {
                int start = bcf->pos - args.flank;
                int end = bcf->pos + args.flank;
                if (start < 0) start = 0;
                int len = 0;
                char *f = faidx_fetch_seq(args.fai,  (char*)hdr->id[BCF_DT_CTG][bcf->rid].key, start, end, &len);
                if (gtf->strand == 1) {
                    char *r = reverse_seq(f, len);
                    free(f);
                    f = r;
                }
                fprintf(args.fasta_fp_out, ">%s:%d-%d|||GN:Z:%s\n%s\n",(char*)hdr->id[BCF_DT_CTG][bcf->rid].key, start, end, dict_name(args.G->gene_name,gtf->gene_name), f);
                free(f);
            }
        }
    }
    if (tmpi) free(tmpi);
    bcf_destroy(bcf);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
    if (out) hts_close(out);
    
    int i, j;
    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j) {
            if (i ==j) continue;
            printf("%c>%c\t%d\t%d\n", "ACGT"[i],"ACGT"[j], geno_stat_ori[i][j], geno_stat[i][j]);
        }
    }

    memory_release();
    
    return 0;
}
