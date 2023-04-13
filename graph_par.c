#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <omp.h>


// matrix size and some others
#include "graph.h"

// some utilities to print matrices, etc ...
#include "util.h"

#define DIFFTEMPS(a,b) (((b).tv_sec - (a).tv_sec) + ((b).tv_usec - (a).tv_usec)/1000000.)


/* 
 * Transitive closure with the Floyd-Warshall algorithm.
 * Works for directed or undirected graphs.
 * https://fr.wikipedia.org/wiki/Algorithme_de_Warshall
 *
 * @param a : the adjacency matrix of the graph
 * @param c : output: the adjacency matrix resulting from tansitive closure   
 */
void warshall(size_t n, int a[n][n], int c[n][n]) {

#ifdef PROGRESS
    size_t percent=0;
#endif
    //---------------------------------------------------------------------------
    // init matrix c from matrix a
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            c[i][j] = a[i][j];

    /*  updates the value of c[i][j] to true if there is a path 
     *  from i to j that goes through k. */
    for (size_t k = 0; k < n; k++) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                c[i][j] = c[i][j] || (c[i][k] && c[k][j]);
            }
        }
#ifdef PROGRESS
        if ((k*100)%n==0) {
            percent++;
            fprintf(stderr, "(%3zu%%)", percent); fflush(stderr);fprintf(stderr, "\b\b\b\b\b\b");
        }
#endif
    }
}

void warshall_par(size_t n, int a[n][n], int c[n][n]) {


     int i, j, k;
 #pragma omp parallel shared(a, c, n) private(i, j, k)
    {

         // init matrix c from matrix a
     #pragma omp for 
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                c[i][j] = a[i][j];
       //* print_matrix(): printed 1000 lines.
      //Init : 0.056649s, Compute : 1.430937s

      // version par
     // * print_matrix(): printed 1000 lines.
    //Init : 0.040848s, Compute : 0.318178s
        
        //  updates the value of c[i][j] to true if there is a path
        //  from i to j that goes through k.
        #pragma omp for
        for (k = 0; k < n; k++) {
            for (i = 0; i < n; i++) {
                for (j = 0; j < n; j++) {
                    c[i][j] = c[i][j] || (c[i][k] && c[k][j]);
                }
            }
        }


         
    }
    

    
}


void warshall_par2(size_t n, int a[n][n], int c[n][n]) {
    
   // split c into 4 parts and init from a
  
    //

    
}

int min(int a, int b) {
    return a < b ? a : b;
}

void blocked_floyd_warshall( int n,int W[][n],int c[][n])
{
    int i, j, k, ii, jj, kk;
    int b = 4;  // block size
    int Wkk[b][b], Wkj[b][b], Wik[b][b], Wij[b][b];
    // init matrix c from matrix W
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            c[i][j] = W[i][j];


    // split W into blocks
    for (kk = 0; kk < n; kk += b) {
        for (ii = 0; ii < n; ii += b) {
            for (jj = 0; jj < n; jj += b) {
                // process block W[kk][kk]
                for (i = ii; i < (ii + b); i++) {
                    for (j = jj; j < (jj + b); j++) {
                        Wkk[i - ii][j - jj] = c[i][j];
                    }
                }
                for (k = kk; (k < kk + b )&&( k < n); k++) {
                    for (i = ii; i < (ii + b) ; i++) {
                        for (j = jj; j < (jj + b); j++) {
                            c[i - ii][j - jj] = Wkk[i - ii][j - jj]|| ( c[i - ii][k - kk] && c[k][j - jj]);
                        }
                    }
                }
                for (i = ii; i < (ii + b ); i++) {
                    for (k = kk; (k < kk + b )&&( k < n); k++) {
                        for (j = jj; j < (jj + b); j++) {
                            Wik[i - ii][j - jj] = c[i - ii][j - jj];
                            Wik[i - ii][j - jj] = Wik[i - ii][j - jj] || (c[i - ii][k - kk] && c[k][j - jj]);
                                                  //c[i][j] || (c[i][k] && c[k][j]);
                        }
                    }
                }
                // for( ; i+o<N*C && i < C*n; i++ )
                //int maxII = min( N*C-o, C*n);
                //int maxII =  min( ii + b, n);
                //ii<maxII 
                // partially dependent phase
                #pragma omp parallel for private(j, k, Wkj)
                for (j = jj; j < (jj + b); j++) {
                    if (j == kk) continue;
                    for (k = kk; (k < kk + b) && (k < n); k++) {
                        for (i = ii; i < (ii + b); i++) {
                            Wkj[i - ii][k - kk] = c[k][j - jj];
                        }
                    }
                    for (i = ii; i < (ii + b); i++) {
                        for (k = kk; (k < kk + b) && (k < n); k++) {
                            c[i - ii][j - jj] = c[i - ii][j - jj]||( Wik[i - ii][k - kk] && Wkj[i - ii][k - kk]);
                        }
                    }
                }
                  
                // partially dependent phase
                #pragma omp parallel for private(i, k, Wik)
                for (i = ii; i <( ii + b); i++) {
                    if (i == kk) continue;
                    for (k = kk;( k < kk + b) && (k < n); k++) {
                        for (j = jj; j < (jj + b); j++) {
                            Wik[j - jj][k - kk] = c[i][k - kk];
                        }
                    }
                    for (j = jj; j < (jj + b); j++) {
                        for (k = kk; (k < kk + b) && (k < n); k++) {
                            c[i - ii][j - jj] = c[i - ii][j - jj] ||(Wik[j - jj][k - kk] && Wkj[i - ii][k - kk]);
                        }
                    }
                }

                # pragma omp parallel for private(i, j, Wij)
                for (i = ii; i < (ii + b); i++) {
                    if (i == kk) continue;
                    for (j = jj; j < (jj + b); j++) {
                        if (j == kk) continue;
                        for (k = kk; (k < kk + b) && (k < n); k++) {
                            Wij[i - ii][j - jj] = c[i][j - jj];
                        }
                    }
                    for (j = jj; j < (jj + b) ; j++) {
                        if (j == kk) continue;
                        for (k = kk; (k < kk + b) &&( k < n); k++) {
                            c[i - ii][j - jj] = c[i - ii][j - jj] ||(Wik[j - jj][k - kk] && Wkj[i - ii][k - kk]);
                        }
                    }
                }
            }
        }
    }
}

/**
 * main
 */
int main (int argc, char **argv) {

	  
	char * filename = 0;
	char * input_type = 0;
	parse_options(argc, argv, &filename, &input_type);

      struct timeval tv_init, tv_begin, tv_end;
      gettimeofday (&tv_init, NULL);


      size_t n = 0;
      // read from an adjacency matrix
      // determine size n of matrix
      if (strcmp(input_type, INPUT_TYPE_ADJ)==0) {
          // first get the adjacency matrix size
          n = matrix_lines_from_file(filename);
          fprintf(stderr, "* %s has %ld lines.\n", filename, n);
      }
      // read from a file, pairs of 2 nodes per line
      // determine size n of matrix
      if (strcmp(input_type, INPUT_TYPE_PAIRS)==0) {
          // first get the max integers in the pairs
          n = find_max(filename, "\t");
          // assume values in file may start at 0, so add one to max value. 
          n++;
          fprintf(stderr, "* max value found in %s: %ld (about to build a %ldx%ld adjacency matrix)\n", filename, n, n, n);

      }
      // malloc a VLA as a pointer to an array of n integers
      int (*a)[n] = malloc(sizeof(int[n][n]));
      if (a==0) {
          fprintf(stderr, "* could not alloc %ld bytes. Abort.\n", n*n);
          exit(1);
      }

      // Now load, from one type or another
      if (strcmp(input_type, INPUT_TYPE_ADJ)==0) {
          // read data from a file containing an adjacency matrix
          matrix_from_adj_file(filename, n, a);
      }
      // read from a file, pairs of 2 nodes per line
      if (strcmp(input_type, INPUT_TYPE_PAIRS)==0) {
          // read data from a file containing a list of pairs
          // : adjust the separator "\t", ",", ... depending on your file
          matrix_from_pairs_file(filename, n, a, "\t");
      }

      int (*c)[n] = malloc(sizeof(int[n][n]));
      if (c==0) {
          fprintf(stderr, "* could not alloc %ld bytes. Abort.\n", n*n);
          exit(1);
      }

#ifdef DEBUG
      print_matrix ("a_orig", n , a, 0, n, 0, n, 0, false);
#endif
      //Compute
       printf("original matrix\n");
        print_matrix (NULL, n, a, 0, n, 0, n, 0, false);
      fprintf(stderr, "* starting computation (n=%ld) ... ", n);
      fflush(stderr);
      gettimeofday (&tv_begin, NULL);
      warshall_par(n, a, c);
      //blocked_floyd_warshall(n, a, c);
      gettimeofday (&tv_end, NULL);
      fprintf(stderr, " done.\n");

      print_matrix (NULL, n, c, 0, n, 0, n, 0, false);

#ifdef CCOMP
      size_t ncomps;
      comp_t * ccomp = make_ccomp_digraph (n, c, &ncomps);

      free(a);
      free(c);
      fprintf (stderr, "* %ld connected components after make_ccomp_digraph.\n", ncomps);

      // write the graph transitivly closed in dot format
      char *output_file = malloc(strlen(OUTPUT_TYPE) + strlen(OUTPUT_EXT) + strlen(filename) + 1); 
      sprintf(output_file, "%s%s%s", filename, OUTPUT_TYPE, OUTPUT_EXT);
      FILE *fgraph_clos = fgraph_clos = fopen (output_file, "w");
      if (fgraph_clos == NULL) {
          fprintf (stderr, "Cannot open file %s for writing.\n", output_file);
          exit (EXIT_FAILURE);
      }
      print_dot (fgraph_clos, ccomp, ncomps);
      fclose (fgraph_clos);
#endif

      //---------------------------------------------------------------------------
      //Execution times
      fprintf (stderr, "Init : %lfs, Compute : %lfs\n",
              DIFFTEMPS (tv_init, tv_begin), DIFFTEMPS (tv_begin, tv_end));

      return 0;
}
