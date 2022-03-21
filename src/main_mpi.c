#if defined(__STDC__)
#  if (__STDC_VERSION__ >= 199901L)
#     define _XOPEN_SOURCE 700
#  endif
#endif
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include <time.h>
#include <kdtree_shared.h>
#include <kdtree_mpi.h>


int main(int argc,char ** argv )
{
  int rank, N_p;
  
  MPI_Init( &argc, &argv );
  MPI_Comm_size(MPI_COMM_WORLD, &N_p);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  if(!IsPowerOfTwo(N_p)) {
    if(rank==0)
      printf("The number of processes (%d) must be a power of 2\n", N_p);
    exit (EXIT_FAILURE);
  }
  
  int N_glo = 0;
  int N = 0;
  struct kd_point_t *points ;
  
#if defined(LOAD_DATA_FROM_FILE)
  load_from_file_distributed_P2P( argc, argv, &N, &N_glo, NDIM, &points, N_p, rank);
#else
  generate_random_points_distributed_P2P( argc, argv, &N, &N_glo, NDIM, &points, N_p, rank);
  
#endif
  
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0) {
    printf("Number of processes=%d\n", N_p);
    printf("Number of points=%d\n", N_glo);
  }
  int axis = -1;
  int level = -1;
#if defined(DEBUG_POINTS)
  show_array_P2P( points, 0, N, level, 'A') ;
#endif
  struct kd_node_t * root;
  
  root = (struct kd_node_t*) malloc (sizeof(struct kd_node_t));
  
  double start, end, maxtime, mintime, avgtime, time;
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();
  build_kdtree_distributed_P2P(root, &points, N, NDIM, axis, level, MPI_COMM_WORLD);
  end = MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);
  time = end-start;
  /*compute max, min, and average timing statistics*/
  MPI_Reduce(&time, &maxtime, 1, MPI_DOUBLE,MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&time, &mintime, 1, MPI_DOUBLE, MPI_MIN, 0,MPI_COMM_WORLD);
  MPI_Reduce(&time, &avgtime, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
  if (rank == 0) {
    avgtime /= N_p;
    printf("[MPI kd tree] Min: %lf Max: %lf Avg: %lf\n", mintime, maxtime,avgtime);
  }
  
#if defined(DEBUG_TREE)
  debug_tree_P2P(root);
#endif
  
#if defined(PRINT_TREE)
  if(rank==0) printf("A) SIZE= %d\n",N_p );
  struct kd_node_t * tmp = root;
#ifndef STORE
  while(tmp->kid) tmp = tmp->kid;
#endif
  if(tmp->split)
    printf("B) %d %3d P%d_%f_%f %f %f\n",  rank, tmp->level, rank,tmp->split->kd_coord[0], tmp->split->kd_coord[1], tmp->split->kd_coord[0], tmp->split->kd_coord[1] );
  int time_to_sleep=1*rank;
  sleep(time_to_sleep);
  print_tree_P2P(root);
#endif
  clear_P2P(root);
  free(points);
  MPI_Finalize();
  return 0;
}
//mpicc MPI_kdtree_A2A_file.c -DDOUBLE_PRECISION -march=native  -O3 -Wall -o MPI_kdtree_A2A_file.x -lm



