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
void ascii_to_bin(const char *asciiFile, const char *binFile) {
    FILE *fp_ascii = fopen(asciiFile, "r");
    FILE *fp_bin = fopen(binFile, "wb");

    if (fp_ascii == NULL) {
        perror(asciiFile);
        exit(1);
    }
  
    if (fp_bin == NULL) {
        perror(binFile);
        exit(1);
    }
  
    struct DumpHeader header;

    // Read and write the header data
    fscanf(fp_ascii, "%la %ld %d %d %d %d", &header.t, &header.len, &header.i, &header.depth, &header.npe, &header.version);

    //fwrite(&header, sizeof(struct DumpHeader), 1, fp_bin);

    // Read and write the MPI dimensions
    fscanf(fp_ascii, "%la %la %la", &header.n.x, &header.n.y, &header.n.z);
    //fwrite(&header.n, sizeof(header.n), 1, fp_bin);
    fwrite(&header, sizeof(struct DumpHeader), 1, fp_bin);

    // Read and write the scalar names
    char name[256];
    fscanf(fp_ascii, "%s", name);

    for (int i = 0; i < header.len; i++) {
        char name[256]; // assuming a maximum length for the scalar names
        fscanf(fp_ascii, "%s", name);
        unsigned len = strlen(name);
        fwrite(&len, sizeof(len), 1, fp_bin);
        fwrite(name, sizeof(char), len, fp_bin);
        fprintf(stdout, "Read scalar name %s %d\n", name, len);
    }

    // Read and write the coordinates
    double o[4];
    fscanf(fp_ascii, "%la %la %la %la", &o[0], &o[1], &o[2], &o[3]);
    fwrite(o, sizeof(double), 4, fp_bin);
  
    unsigned flags;
    double val;

    // Notice the change in the loop condition here
    while (fscanf(fp_ascii, "%u", &flags) == 1) {
        fwrite(&flags, sizeof(flags), 1, fp_bin);
      
        for (int i = 0; i < header.len; i++) {
            fscanf(fp_ascii, "%la", &val);
            fwrite(&val, sizeof(val), 1, fp_bin);
        }
    }

    fclose(fp_ascii);
    fclose(fp_bin);
}
#else
void ascii_to_bin(const char *asciiFile, const char *binFile) {
    FILE *fp_ascii = fopen(asciiFile, "r");
    FILE *fp_bin = fopen(binFile, "wb");

    if (fp_ascii == NULL) {
        perror(asciiFile);
        exit(1);
    }
  
    if (fp_bin == NULL) {
        perror(binFile);
        exit(1);
    }
  
    struct DumpHeader header;

    // Read and write the header data
    fscanf(fp_ascii, "%le %ld %d %d %d %d", &header.t, &header.len, &header.i, &header.depth, &header.npe, &header.version);

    //fwrite(&header, sizeof(struct DumpHeader), 1, fp_bin);

    // Read and write the MPI dimensions
    fscanf(fp_ascii, "%le %le %le", &header.n.x, &header.n.y, &header.n.z);
    //fwrite(&header.n, sizeof(header.n), 1, fp_bin);
    fwrite(&header, sizeof(struct DumpHeader), 1, fp_bin);

    // Read and write the scalar names
    char name[256];
    fscanf(fp_ascii, "%s", name);

    for (int i = 0; i < header.len; i++) {
        char name[256]; // assuming a maximum length for the scalar names
        fscanf(fp_ascii, "%s", name);
        unsigned len = strlen(name);
        fwrite(&len, sizeof(len), 1, fp_bin);
        fwrite(name, sizeof(char), len, fp_bin);
        fprintf(stdout, "Read scalar name %s %d\n", name, len);
    }

    // Read and write the coordinates
    double o[4];
    fscanf(fp_ascii, "%le %le %le %le", &o[0], &o[1], &o[2], &o[3]);
    fwrite(o, sizeof(double), 4, fp_bin);
  
    unsigned flags;
    double val;

    // Notice the change in the loop condition here
    while (fscanf(fp_ascii, "%u", &flags) == 1) {
        fwrite(&flags, sizeof(flags), 1, fp_bin);
      
        for (int i = 0; i < header.len; i++) {
            fscanf(fp_ascii, "%le", &val);
            fwrite(&val, sizeof(val), 1, fp_bin);
        }
    }

    fclose(fp_ascii);
    fclose(fp_bin);
}
#endif

int main(int argc, char * argv[])
{
    if(argc != 3)
        fprintf(stderr, "2 arguments are required.\n");

    ascii_to_bin(argv[1], argv[2]);

    return 0;
}