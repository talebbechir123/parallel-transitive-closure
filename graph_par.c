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

    // set the number of threads
    //omp_set_num_threads(6);

 #pragma omp parallel shared(a, c, n) private(i, j, k) 
    {

        #pragma omp for schedule(static)
        for (k = 0; k < n; k++) {
            for (i = 0; i < n; i++) {
                for (j = 0; j < n; j++) {
                    c[i][j] = c[i][j] || (c[i][k] && c[k][j]);
                }
            }
        }


         
    }
    

    
}


void warshall_col(size_t n, int a[n][n], int c[n][n]){
  
    int i, j, k;
    // init matrix c from matrix a
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            c[i][j] = a[i][j];
    
    // updates the value of c[i][j] to true if there is a path
    #pragma omp parallel  shared(c, n)  
    {
       
        #pragma omp for schedule(guided,50) private(i, j)
        for (size_t k = 0; k < n; k++) {
            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    c[i][j] = c[i][j] || (c[i][k] && c[k][j]);
                }
            }
        }

    }

}

//use MPI to parallelize the warshall algorithm


int min(int a, int b) {
    return a < b ? a : b;
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
      //use MPI to parallelize the warshall algorithm
      int numprocs, myid;
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        MPI_Status status;
        // determine the number of rows each process will handle
        int rows_per_proc = n / numprocs;
        // determine the number of rows that will be handled by the last process
        int last_rows = n - (numprocs - 1) * rows_per_proc;
        // determine the starting row for each process
        int start_row = myid * rows_per_proc;
        // determine the ending row for each process
        int end_row = start_row + rows_per_proc - 1;
        // determine the starting row for the last process
        //let process 0 initialize the matrix
       

      printf("original matrix\n");
      print_matrix (NULL, n, a, 0, n, 0, n, 0, false);
      fprintf(stderr, "* starting computation (n=%ld) ... ", n);
      fflush(stderr);
      gettimeofday (&tv_begin, NULL);
     //warshall_par(n, a, c);
     // warshall_col(n, a, c);

       if (myid == 0) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    c[i][j] = a[i][j];
                }
            }
        }

        // broadcast the matrix to all processes after process 0
        if (myid != 0) {
            MPI_Recv(&c, n * n, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        }
        else {
            for (int i = 1; i < numprocs; i++) {
                MPI_Send(&c, n * n, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
        }
        // determine the starting row for each process
        int start_row = myid * rows_per_proc;
        // determine the ending row for each process
        int end_row = start_row + rows_per_proc - 1;
        // each process will compute the shortest path for the rows it is responsible for
        for (int k = 0; k < n; k++) {
            for (int i = start_row; i <= end_row; i++) {
                for (int j = 0; j < n; j++) {
                    c[i][j] = c[i][j] || (c[i][k] && c[k][j]);
                }
            }
            // process 0 will send the updated rows to the other processes
            if (myid == 0) {
                for (int i = 1; i < numprocs; i++) {
                    MPI_Send(&c[i * rows_per_proc], rows_per_proc * n, MPI_INT, i, 0, MPI_COMM_WORLD);
                }
            }
            // all other processes will receive the updated rows from process 0
            else {
                MPI_Recv(&c[start_row], rows_per_proc * n, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            }
        }

        // process 0 will gather the results from all the other processes
        if (myid == 0) {
            for (int i = 1; i < numprocs; i++) {
                MPI_Recv(&c[i * rows_per_proc], rows_per_proc * n, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
            }
        }

        // all other processes will send their results to process 0
        else {
            MPI_Send(&c[start_row], rows_per_proc * n, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }


       

       
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
