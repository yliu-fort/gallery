/**
# Convert dump file to binary/ascii format

Convert dump file produced by dump() in Basilisk to binary/ascii format.
Use %.12le format for human readable floating point expression or
use %la format for exact floating point conversion so that the reconstructed
binary file is bit-wise correct.
Author: Yuxuan Liu
Date: 2023/05/21
*/

#include "common.h"

struct DumpHeader {
  double t;
  long len;
  int i, depth, npe, version;
  coord n;
};

#ifdef EXACT_FP
void bin_to_ascii(const char *binFile, const char *asciiFile) {
    FILE *fp_bin = fopen(binFile, "rb");
    FILE *fp_ascii = fopen(asciiFile, "w");

    if (fp_bin == NULL) {
        perror(binFile);
        exit(1);
    }
  
    if (fp_ascii == NULL) {
        perror(asciiFile);
        exit(1);
    }
  
    struct DumpHeader header;
    fread(&header, sizeof(header), 1, fp_bin);

    // Read and write the header data
    fprintf(fp_ascii, "%la %ld %d %d %d %d\n", header.t, header.len, header.i, header.depth, header.npe, header.version);

    // Read and write the MPI dimensions
    fprintf(fp_ascii, "%la %la %la\n", header.n.x, header.n.y, header.n.z);

    // Read and write the scalar names
    fprintf(fp_ascii, "isleaf ");

    for (int i = 0; i < header.len; i++) {
        unsigned len;
        fread(&len, sizeof(len), 1, fp_bin);
        char *name = malloc((len + 1) * sizeof(char));
        fread(name, sizeof(char), len, fp_bin);
        name[len] = '\0';
        fprintf(fp_ascii, "%s ", name);
        free(name);
    }

    // Read and write the coordinates
    double o[4];
    fread(o, sizeof(double), 4, fp_bin);
    fprintf(fp_ascii, "%la %la %la %la\n", o[0], o[1], o[2], o[3]);
  
    unsigned flags;
    int cell_size = sizeof(unsigned) + header.len*sizeof(double);
    long pos = ftell(fp_bin);
    int n;
  
    while ((n = fread(&flags, sizeof(flags), 1, fp_bin)) == 1) {
        fseek(fp_bin, -sizeof(unsigned), SEEK_CUR);
        unsigned flags;
        fread(&flags, sizeof(flags), 1, fp_bin);
        fprintf(fp_ascii, "%u ", flags);
      
        double val;
        for (int i = 0; i < header.len; i++) {
            fread(&val, sizeof(val), 1, fp_bin);
            fprintf(fp_ascii, "%la ", val);
        }
        fprintf(fp_ascii, "\n");

        pos += cell_size;
        fseek(fp_bin, pos, SEEK_SET);
    }

    fclose(fp_bin);
    fclose(fp_ascii);
}
#else
void bin_to_ascii(const char *binFile, const char *asciiFile) {
    FILE *fp_bin = fopen(binFile, "rb");
    FILE *fp_ascii = fopen(asciiFile, "w");

    if (fp_bin == NULL) {
        perror(binFile);
        exit(1);
    }
  
    if (fp_ascii == NULL) {
        perror(asciiFile);
        exit(1);
    }
  
    struct DumpHeader header;
    fread(&header, sizeof(header), 1, fp_bin);

    // Read and write the header data
    fprintf(fp_ascii, "%.12le %ld %d %d %d %d\n", header.t, header.len, header.i, header.depth, header.npe, header.version);

    // Read and write the MPI dimensions
    fprintf(fp_ascii, "%.12le %.12le %.12le\n", header.n.x, header.n.y, header.n.z);

    // Read and write the scalar names
    fprintf(fp_ascii, "isleaf ");

    for (int i = 0; i < header.len; i++) {
        unsigned len;
        fread(&len, sizeof(len), 1, fp_bin);
        char *name = malloc((len + 1) * sizeof(char));
        fread(name, sizeof(char), len, fp_bin);
        name[len] = '\0';
        fprintf(fp_ascii, "%s ", name);
        free(name);
    }

    // Read and write the coordinates
    double o[4];
    fread(o, sizeof(double), 4, fp_bin);
    fprintf(fp_ascii, "%.12le %.12le %.12le %.12le\n", o[0], o[1], o[2], o[3]);
  
    unsigned flags;
    int cell_size = sizeof(unsigned) + header.len*sizeof(double);
    long pos = ftell(fp_bin);
    int n;
  
    while ((n = fread(&flags, sizeof(flags), 1, fp_bin)) == 1) {
        fseek(fp_bin, -sizeof(unsigned), SEEK_CUR);
        unsigned flags;
        fread(&flags, sizeof(flags), 1, fp_bin);
        fprintf(fp_ascii, "%u ", flags);
      
        double val;
        for (int i = 0; i < header.len; i++) {
            fread(&val, sizeof(val), 1, fp_bin);
            fprintf(fp_ascii, "%.12le ", val);
        }
        fprintf(fp_ascii, "\n");

        pos += cell_size;
        fseek(fp_bin, pos, SEEK_SET);
    }

    fclose(fp_bin);
    fclose(fp_ascii);
}
#endif

int main(int argc, char * argv[])
{
    if(argc != 3)
        fprintf(stderr, "2 arguments are required.\n");

    bin_to_ascii(argv[1], argv[2]);

    return 0;
}