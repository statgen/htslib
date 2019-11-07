
#include "../../htslib/vcf.h"

int main(int argc, char** argv)
{
  int exit_code = 0;
  int use_sav = 1; //argc > 1 && argv[1][0] != '0' ? 1 : 0;

  htsFile* in = hts_open("../test_file.vcf", "r");
  //htsFile* in = hts_open("/net/topmed8/working/lefaivej/freeze7/phasing/ligated-bcf-chunks/topmed.chr10.gtonly.filtered.rehdr.pass_only.phased.ligated.100000001-101000000.bcf", "r");
  htsFile* out = hts_open(use_sav ? "../test_file_out.sav" : "../test_file_out.bcf", "wb");

  bcf_hdr_t* hdr = bcf_hdr_read(in);
  bcf1_t* rec = bcf_init1();

  if (!hdr)
  {
    fprintf(stderr, "Error reading header\n");
    exit_code = 1;
  }

  if (!rec)
  {
    fprintf(stderr, "Error initializing\n");
    exit_code = 1;
  }

//  bcf_hdr_append(hdr, "##INFO=<ID=GTV,Number=.,Type=Integer,Description=\"GT Values\">");
//  bcf_hdr_append(hdr, "##INFO=<ID=GTO,Number=.,Type=Integer,Description=\"GT Offsets\">");

  if (exit_code == 0)
  {
    if (bcf_hdr_write(out, hdr) < 0)
    {
      fprintf(stderr, "Error writing header.\n");
      exit_code = 1;
    }
    else
    {

      int* sbuf_val = malloc(sizeof(int) * bcf_hdr_nsamples(hdr) * 2);
      unsigned int* sbuf_off = malloc(sizeof(int) * bcf_hdr_nsamples(hdr) * 2);
      int* buf = NULL;
      int buf_n = 0;

      int res = 0;
      while ((res = bcf_read1(in, hdr, rec)) >= 0)
      {
        bcf_unpack(rec, BCF_UN_FMT);
        bcf_update_format(hdr, rec, "GQ", NULL, 0, 0);
        bcf_update_format(hdr, rec, "DP", NULL, 0, 0);
        bcf_update_format(hdr, rec, "HQ", NULL, 0, 0);
        bcf_update_format(hdr, rec, "HDS", NULL, 0, 0);

        bcf_get_format_int32(hdr, rec, "GT", &buf, &buf_n);
        int sbuf_n = 0;
        for (int i = 0; i < buf_n; ++i)
        {
          if (bcf_gt_allele(buf[i]))
          {
            sbuf_val[sbuf_n] = bcf_gt_allele(buf[i]);
            sbuf_off[sbuf_n++] = i;
          }
        }

        if (use_sav)
          sav_update_format(hdr, rec, "GT", sbuf_val, buf_n, BCF_HT_INT, sbuf_off, sbuf_n);
        else
          bcf_update_genotypes(hdr, rec, NULL, 0); //buf, buf_n);

//        bcf_update_info_int32(hdr, rec, "GTV", sbuf_val, sbuf_n);
//        bcf_update_info_int32(hdr, rec, "GTO", sbuf_off, sbuf_n);

        if (bcf_write1(out, hdr, rec) < 0)
        {
          fprintf(stderr, "Error writing output.\n");
          exit_code = 1;
          break;
        }
      }

      if (res < -1)
      {
        fprintf(stderr, "Error reading input.\n");
        exit_code = 1;
      }
    }
  }

  //bcf_close(in);
  //bcf_close(out);

  if (rec)
    bcf_destroy1(rec);
  if (hdr)
    bcf_hdr_destroy(hdr);

  bcf_close(in);
  bcf_close(out);

  return exit_code;
}