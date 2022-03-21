#if defined(__STDC__)
#  if (__STDC_VERSION__ >= 199901L)
#     define _XOPEN_SOURCE 700
#  endif
#endif

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <kdtree_shared.h>

#if defined(_OPENMP)
#define CPU_TIME (clock_gettime( CLOCK_REALTIME, &ts ), (double)ts.tv_sec + \
                  (double)ts.tv_nsec * 1e-9)

#define CPU_TIME_th (clock_gettime( CLOCK_THREAD_CPUTIME_ID, &myts ), (double)myts.tv_sec + \
                     (double)myts.tv_nsec * 1e-9)

#else

#define CPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec + \
                  (double)ts.tv_nsec * 1e-9)
#endif

int main(int argc,char ** argv )
{
  int N = 1;
  struct kd_point_t *points;
  
#if defined(LOAD_DATA_FROM_FILE)
  load_from_file_shared( argc, argv, &N, NDIM, & points);
#else
  generate_random_points_shared(argc, argv, &N, NDIM, & points);
#endif
  
  printf("Number of points=%d\n", N);
  struct kd_node_t *root;
  root = (struct kd_node_t*) malloc (sizeof(struct kd_node_t));
  int level = -1;
  int axis = -1;
#if defined(DEBUG_POINTS)
  show_array_shared( points, 0, N, level, 'A') ;
#endif
  struct  timespec ts;
  double wtstart = CPU_TIME;
  build_kdtree_serial(root, points, N, NDIM, axis, level);
  double wtend = CPU_TIME;
  printf("Serial kd tree: %g sec elapsed\n", wtend - wtstart);
#if defined(DEBUG_TREE)
  debug_tree_shared(root);
#endif
  
#if defined(PRINT_TREE)
  printf("A) SIZE= 1\n" );
  struct kd_node_t * tmp = root;
  if(tmp->split)
    printf("B) 0 %3d  P0_%f_%f %f %f\n", tmp->level, tmp->split->kd_coord[0], tmp->split->kd_coord[1],tmp->split->kd_coord[0], tmp->split->kd_coord[1] );
  print_tree_shared(root);
#endif
  clear_shared(root);
  free(points);
  return 0;
}


// gcc-11 Serial_kdtree_final.c -O3 -Wall -o Serial_kdtree_final.x
// gcc Serial_kdtree_final.c -O3 -Wall -o Serial_kdtree_final.x

// ./Serial_kdtree_final.x 10.input
// leaks -atExit -- ./Serial_kdtree_final.x 10.input
