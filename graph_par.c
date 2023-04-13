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
     for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                c[i][j] = a[i][j];

 #pragma omp parallel shared(a, c, n) private(i, j, k)
    {

         // init matrix c from matrix a
     
        
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



int min(int a, int b) {
    return a < b ? a : b;
}

void floyd_w2(int b,int c[b][b],int a[b][b],int B[b][b])
{
    int i, j, k;
    for (k = 0; k < b; k++) {
        for (i = 0; i < b; i++) {
            for (j = 0; j < b; j++) {
                c[i][j] = min(c[i][j], a[i][k] + B[k][j]);
                 
            }
        }
    }
    
}

    
// Blocked Floyd-Warshall algorithm with OpenMP
void blocked_floyd_warshall(int W[][n], int n) 
{
    int B = 4;  // number of blocks
    int k, i, j;

    // Loop over blocks
    #pragma omp parallel for schedule(static)
    for (k = 0; k < B; k++) {
        int Wkk[b][b], Wkj[b][b], Wik[b][b], Wij[b][b];
        int W_kk[b][b], W_kj[b][b], W_ik[b][b], W_ij[b][b];

        // Copy the block W_kk into a local buffer
        for (i = 0; i < b; i++) {
            for (j = 0; j < b; j++) {
                Wkk[i][j] = W[k*b+i][k*b+j];
            }
        }

        // Dependent phase
        floyd(Wkk, Wkk, Wkk, b);

        // Partially dependent phase
        #pragma omp for schedule(static)
        for (j = 0; j < B; j++) {
            if (j != k) {
                // Copy the block W_kj into a local buffer
                for (i = 0; i < b; i++) {
                    for (int l = 0; l < b; l++) {
                        Wkj[i][l] = W[k*b+i][j*b+l];
                    }
                }
                floyd(Wkj, Wkk, Wkj, b);
                // Copy the local buffer back into W_kj
                for (i = 0; i < b; i++) {
                    for (int l = 0; l < b; l++) {
                        W[k*b+i][j*b+l] = Wkj[i][l];
                    }
                }
            }
        }

        #pragma omp for schedule(static)
        for (i = 0; i < B; i++) {
            if (i != k) {
                // Copy the block W_ik into a local buffer
                for (int l = 0; l < b; l++) {
                    for (j = 0; j < b; j++) {
                        Wik[l][j] = W[i*b+l][k*b+j];
                    }
                }
                floyd(Wik, Wik, Wkk, b);
                // Copy the local buffer back into W_ik
                for (int l = 0; l < b; l++) {
                    for (j = 0; j < b; j++) {
                        W[i*b+l][k*b+j] = Wik[l][j];
                    }
                }
            }
        }

        // Independent phase
        #pragma omp for schedule(static)
        for (i = 0; i < B; i++) {
            if (i != k) {
                // Copy the block W_ij into a local buffer 

            for (int l = 0; l < b; l++) {
                for (j = 0; j < B; j++) {
                    if (j != k) {
                        Wij[l][j] = W[i*b+l][j*b+k];
                    }
                }
            }
            floyd(Wij, Wik, Wkj, b);
            // Copy the local buffer back into W_ij
            for (int l = 0; l < b; l++) {
                for (j = 0; j < B; j++) {
                    if (j != k) {
                        W[i*b+l][j*b+k] = Wij[l][j];
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
     // warshall_par(n, a, c);
        blocked_floyd_warshall(n, a, c);
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
