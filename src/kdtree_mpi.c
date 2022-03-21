#if defined(__STDC__)
#  if (__STDC_VERSION__ >= 199901L)
#     define _XOPEN_SOURCE 700
#  endif
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mpi.h"
#include <stdbool.h>
#include "kdtree_shared.h"
#include "kdtree_mpi.h"

bool IsPowerOfTwo(int x)
{
  return (x != 0) && ((x & (x - 1)) == 0);
}

void glb_mean(// input
	      const kd_point *points,
	      const int numof_points,
	      const int glb_numof_points,
	      const int axis,
	      MPI_Comm comm,
	      // output
	      float_t *glb_mu)
/*
 * glb_mean returns the global mean (glb_mu) of the points in the MPI Communicator comm
 * according to the coordinate specified by axis.
 * numof_points are the number of points that belong to the process
 * glb_numof_points are the number of points that belong to all tha process in comm
 */
{
  float_t mu = 0.0;
  for(int i = 0; i < numof_points; i++) {
    mu += points[i].kd_coord[axis];
  }
  MPI_Allreduce(&mu, glb_mu, 1, MPI_float_t, MPI_SUM, comm);
  *glb_mu /= (float_t)glb_numof_points;
}


void find_median(// input
		 const kd_point * arr,
		 const int N,
		 const int ks,
		 const int axis,
		 MPI_Comm comm,
		 // output
		 float_t * median)
/*
 * find_median returns find the ks-th smallest element,
 * according to the direction specified by axis,
 * among all the points in the MPI Communicator comm.
 * Each process passes to find_median its own N points in the array of points arr.
 * To get the median call find_median with ks = glb_N/2
 * where glb_N is the total number of points in comm
 * Basically a parallel quick select
 */
{
  kd_point * S_less;
  S_less =(kd_point*) malloc(N*sizeof(kd_point));
  kd_point * S_great;
  S_great =(kd_point*) malloc(N*sizeof(kd_point));
  int glb_N;
  MPI_Allreduce(&N, &glb_N, 1, MPI_INT, MPI_SUM, comm);
  glb_mean(arr, N, glb_N, axis, comm, median);
  int N_less = 0;
  int N_great = 0;
  for(int i = 0; i < N; i++) {
    if(arr[i].kd_coord[axis] > *median)
      S_great[N_great++]=arr[i];
    else
      S_less[N_less++]= arr[i];
  }
  int glb_N_less, glb_N_great;
  MPI_Allreduce(&N_less, &glb_N_less, 1, MPI_INT, MPI_SUM, comm);
  MPI_Allreduce(&N_great, &glb_N_great, 1, MPI_INT, MPI_SUM, comm);
  if( glb_N_less == ks || glb_N == 1 || glb_N == glb_N_less || glb_N == glb_N_great ) {
    free(S_less);
    free(S_great);
    return;
  }
  else if(glb_N_less > ks) {
    find_median( S_less, N_less, ks, axis, comm, median);
    free(S_less);
    free(S_great);
    return ;
  }
  else {
    find_median( S_great, N_great, ks-glb_N_less, axis, comm, median );
    free(S_less);
    free(S_great);
    return;
  }
}


void show_array_P2P(const  kd_point *data,
                    const int start,
                    const int end,
                    const int level,
                    const char c)
/*
 * Function to print the points in data from start to end,
 * according to the level in the tree
 */
{
  int global_size, global_rank;
  MPI_Comm_size(MPI_COMM_WORLD, & global_size);
  MPI_Comm_rank(MPI_COMM_WORLD, & global_rank);
  for ( int i = start; i < end; i++ )
    //printf("p:[%d](%c){%d}-(%d/%d) [%.3f, %.3f],\n", level, c, global_rank, i, end, data[i].kd_coord[0], data[i].kd_coord[1]);
    printf("p:[%d](%c){%d} [%.3f, %.3f],\n", level, c, global_rank, data[i].kd_coord[0], data[i].kd_coord[1]);
}


void Pointwise_exchange(// input:
                        kd_point * X_l, const int N_l,
                        kd_point * X_r, const  int N_r,
                        int ndim,
                        MPI_Comm comm,
                        // output:
                        kd_point ** X_n,
                        int * N_n)
/*
 * Given a left array of points N_l with N_l points and
 * Given a right array of points N_r with N_r points
 * Given a pointer X_n to an array to be resized and modified
 * this function finds a companion patner and performs points exchange
 * modifying the X_n assigning to it N_n points
 */
{
  
  int rank, size;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);
  // Package the points in MPI type to be sent around
  MPI_Datatype MPI_kd_point, tmp_type;
  int array_of_blocklengths[] = {ndim};
  MPI_Aint array_of_displacements[] = {0};
  MPI_Datatype array_of_types[] = {MPI_float_t};
  MPI_Aint lb, extent;
  MPI_Type_create_struct(1, array_of_blocklengths, array_of_displacements, array_of_types, &tmp_type);
  MPI_Type_get_extent(tmp_type, &lb, &extent);
  MPI_Type_create_resized(tmp_type, lb, extent, &MPI_kd_point);
  MPI_Type_commit(&MPI_kd_point);
  //
  MPI_Status status;
  // Size left and right points arrays of the partner
  int N_l_partner;
  int N_r_partner;
  // Some tag, just to be sure the sender and receiver match
  int tag1 = 10;
  int tag2 = 20;
  int tag3 = 30;
  int tag4 = 40;
  // Here, we go and let each other know their dimensions
  if( rank < size/2) { //Process belongs to left child
    MPI_Ssend(&N_l, 1, MPI_INT, size - rank - 1, tag1, comm);
    MPI_Recv(&N_l_partner, 1, MPI_INT, size - rank - 1 , tag2, comm, &status);
    MPI_Ssend(&N_r, 1, MPI_INT, size - rank - 1, tag3, comm);
    MPI_Recv(&N_r_partner, 1, MPI_INT, size - rank - 1, tag4, comm, &status);
  }
  else {
    MPI_Recv(&N_l_partner, 1, MPI_INT, size - rank - 1, tag1, comm, &status);
    MPI_Ssend(&N_l, 1, MPI_INT, size - rank - 1, tag2, comm);
    MPI_Recv(&N_r_partner, 1, MPI_INT, size - rank - 1 , tag3, comm, &status);
    MPI_Ssend(&N_r, 1, MPI_INT, size - rank - 1, tag4, comm);
  }
  // And here we exchange the points array and create the new one
  if( rank < size/2 ) { //Process belongs to left child
    *N_n = N_l+N_l_partner;
    *X_n =  realloc(*X_n,(*N_n )*sizeof(kd_point));
    MPI_Ssend(X_r, N_r, MPI_kd_point, size - rank - 1, tag1, comm);
    MPI_Recv(*X_n, N_l_partner, MPI_kd_point, size - rank - 1, tag2, comm, &status);
    for (int i = 0; i< N_l; ++i)
      (*X_n)[i+N_l_partner] = X_l[i];
  }
  else {
    *N_n = N_r+N_r_partner;
    *X_n = realloc(*X_n,(*N_n )*sizeof(kd_point));
    MPI_Recv(*X_n, N_r_partner, MPI_kd_point, size - rank - 1, tag1, comm, &status);
    MPI_Ssend(X_l, N_l, MPI_kd_point, size - rank - 1, tag2, comm);
    for (int i = 0; i< N_r; ++i)
      (*X_n)[i+N_r_partner] = X_r[i];
  }
  // Let's free the MPI type
  MPI_Type_free(&MPI_kd_point);
  MPI_Type_free(&tmp_type);
  // and return
  return;
}

void build_kdtree_distributed_P2P(struct kd_node_t * node, struct kd_point_t **points_ptr, const int N, const int ndim, const int axis,const int level, MPI_Comm comm)
{
  // local rank, size in local communicator
  // global rank, global in MPI_COMM_WORLD
  int rank, N_p, world_rank, world_size;
  
  MPI_Comm_size(comm, &N_p);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  if (N<1) {
    fprintf(stderr, "%s", "Error N<1!\n");
    MPI_Abort(MPI_COMM_SELF, MPI_ERR_OTHER);
    return;
  }
  // axis is the splitting dimension from the previous call
  // level is the level of depth from the previous call
  if ( level == MAX_MPI_LEVEL || N_p < 2 || N < MIN_NUM_OF_POINTS ) {
    node->kid = NULL;
#if defined(DEBUG_POINTS)
    show_array_P2P( *points_ptr, 0, N, level, 'Z') ;
#endif
    build_kdtree_serial(node, *points_ptr,  N,  ndim,  axis, level);
    return;
  }
  int myaxis = choose_splitting_dimension ( ndim, axis);
  int mylevel = level + 1;
  node->axis = myaxis;
  node->level = mylevel;
  // Number of points in the communicator
  int glb_N;
  MPI_Allreduce(&N, &glb_N, 1, MPI_INT, MPI_SUM, comm);
  // Global median
  float_t median;
  find_median( *points_ptr, N, glb_N/2,  myaxis, comm, &median);
  
  // Splitting points
  struct kd_point_t * splitting_point;
  splitting_point = NULL;
  struct kd_point_t * splitting_point_left;
  splitting_point_left = NULL;
  struct kd_point_t * splitting_point_right;
  splitting_point_right = NULL;
  int splitting_index = 0;
  int splitting_index_left = 0;
  int splitting_index_right = 0;
  // Dimension of the left and right points array
  int N_l=0;
  int N_r=0;
  // We need to store N points just to be sure
  struct kd_point_t *left_points = (struct kd_point_t*) malloc(N*sizeof(struct kd_point_t));
  struct kd_point_t *right_points = (struct kd_point_t*) malloc(N*sizeof(struct kd_point_t));
  if (!left_points || !right_points)
    {
      fprintf(stderr, "Unable to allocate left and right space\n");
      MPI_Abort(MPI_COMM_SELF, MPI_ERR_NO_MEM);
      return;
    }
  // Useful to fing the splitting points
  float_t tmp_left = 0;
  float_t tmp_right = 0;
  // Lets scan the points to fill the left and right array,
  // and find the splitting points as the maximum of the left array
  // or the minimum of the right array
  for(int n = 0; n < N; ++n) {
    if((*points_ptr)[n].kd_coord[myaxis] <= median ) {
      left_points[N_l]= (*points_ptr)[n];
      if((*points_ptr)[n].kd_coord[myaxis] > tmp_left ) {
	tmp_left = (*points_ptr)[n].kd_coord[myaxis] ;
	splitting_point_left = *points_ptr + n ;
	splitting_index_left = N_l;
      }
      N_l++;
    }
    else {
      right_points[N_r] = (*points_ptr)[n];
      if((*points_ptr)[n].kd_coord[myaxis] < tmp_right ) {
	tmp_right = (*points_ptr)[n].kd_coord[myaxis] ;
	splitting_point_right = *points_ptr + n ;
	splitting_index_right = N_r;
      }
      N_r++;
    }
  }
  // If there are points in the left array,
  // the splitting point is associated to the maximum of the left array
  // We also remove the splitting point from the array of points
  if(N_l>0) {
    // Here we swap the splitting point to be the last element of the left array
    // If you wish to keep the splitting point in the parallel phase comment this line (*)
#if defined (STORE)
    swap(&left_points[N_l-1], &left_points[splitting_index_left] );
    // Here we ``remove'' the last element, i.e. the splitting point
    // and this line (*)
    --N_l;
#endif
    // Here we identify the true splitting point with the left wplitting point
    splitting_point=splitting_point_left;
    // and we save its index
    splitting_index=splitting_index_left;
  }
  else { // In case N_l==0, we do the same for the right array
    // If you wish comment this line (*)
#if defined (STORE)
    swap(&right_points[N_r-1], &right_points[splitting_index_right] );
    // and this line (*)
    --N_r;
#endif
    splitting_point=splitting_point_right;
    splitting_index=splitting_index_right;
  }
  // Honestly, I find myself more confortable having at least one element in each array
  // So in case, either the left or the right array is empty,
  // we just move one point from one array to the other.
  // I know it's maybe not needed at all, I still have to perform check..
  // For sure it's not the right rhing to do, but it shall be a very patological case
  // But let's do this for the time being and we'll check when we'll have time
  /*
    if(N_l==0 && N_r>0) {
    left_points[0] = right_points[N_r-1];
    --N_r;
    ++N_l;
    }
    if(N_r==0&& N_l>0) {
    right_points[0] = left_points[N_l-1];
    --N_l;
    N_r++;
    }
  */
  /*
    printf("[level=%d]-(local_size=%d,local_rank=%d)-(world_size=%d,world_rank=%d):{N=%d:[N_l=%d,N_r=%d]} (* axis=%d *) -> median: (%lf, %lf, %lf)\n",
    mylevel, N_p, rank, world_size, world_rank, N, N_l, N_r, myaxis,
    median, splitting_point->kd_coord[0], splitting_point->kd_coord[1]);
  */
  // If everything is correct and splitting point exists..
  if(splitting_point) {
    // Here, just to be sure, since we reallocate the points array, free and allocate memory,
    // better not having a pointer to the splitting point in the array,
    // but rather we save the splitting point in the node,
    // we'll have to remember to free it later
    // Also, we may or may not store the points in the parallel section
    node->split = (kd_point*) malloc(sizeof(kd_point));
    // If you want to store the data only in the serial part,
    // and just move the points in the parallel section comment this line
#if defined (STORE)
    *node->split = *splitting_point;
#else
    node->split = NULL;
#endif
    // and use this one instead
    
    // new size of the point array
    int new_N;
    // Finally, we Pointwise_exchange
    Pointwise_exchange( left_points, N_l, right_points, N_r, ndim, comm, points_ptr,&new_N );
    // and free the left and right points
    free(left_points);
    free(right_points);
    left_points=NULL;
    right_points=NULL;
    // We are going to split the communicator
    int new_rank, new_size, color;
    MPI_Comm  new_comm;
    color = (rank < N_p/2) ? 0 : 1;
    MPI_Comm_split(comm, color, rank, &new_comm);
    MPI_Comm_size(new_comm, &new_size);
    MPI_Comm_rank(new_comm, &new_rank);
    // We are in the parallel part, we have only one kid
    node->left = NULL;
    node->right = NULL;
    // We allocate the memory for the kid here
    node->kid =(struct kd_node_t*) malloc (sizeof(struct kd_node_t));
    // We recursively build the distributed tree
    build_kdtree_distributed_P2P(node->kid, points_ptr, new_N, ndim, myaxis, mylevel, new_comm);
  }
  if (comm!=MPI_COMM_WORLD)
    MPI_Comm_free(&comm);
  return ;
}



void print_tree_P2P(const struct kd_node_t *node)
/*
 * Function to visualize the tree generated by node, mainly for debugging/checking
 */
{
  if(node==NULL) {return;}
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(node->kid==NULL && node->left==NULL && node->right==NULL) {
    printf("C) %d %3d P%d_%f_%f %f %f\n", rank,node->level, rank, node->split->kd_coord[0], node->split->kd_coord[1] ,node->split->kd_coord[0], node->split->kd_coord[1]);
    return;
  }
  if(node->kid) {
    if(node->split && node->kid->split) {
      printf("C) %d %3d P%d_%f_%f %f %f\n", rank,node->level, rank,node->split->kd_coord[0], node->split->kd_coord[1],node->split->kd_coord[0], node->split->kd_coord[1] );
      printf("D) %d %3d P%d_%f_%f P%d_%f_%f\n", rank,node->level, rank, node->split->kd_coord[0], node->split->kd_coord[1] , rank , node->kid->split->kd_coord[0], node->kid->split->kd_coord[1]);
    }
    print_tree_P2P(node->kid);
  }
  else {
    if (node->left) {
      printf("C) %d %3d P%d_%f_%f %f %f\n", rank,node->level, rank,node->split->kd_coord[0], node->split->kd_coord[1],node->split->kd_coord[0], node->split->kd_coord[1] );
      printf("D) %d %3d P%d_%f_%f P%d_%f_%f\n", rank,node->level, rank, node->split->kd_coord[0], node->split->kd_coord[1] , rank , node->left->split->kd_coord[0], node->left->split->kd_coord[1]);
      print_tree_P2P(node->left);
    }
    
    if (node->right) {
      printf("C) %d %3d P%d_%f_%f %f %f\n", rank, node->level, rank, node->split->kd_coord[0], node->split->kd_coord[1],node->split->kd_coord[0], node->split->kd_coord[1] );
      printf("D) %d %3d P%d_%f_%f P%d_%f_%f\n", rank,node->level, rank, node->split->kd_coord[0], node->split->kd_coord[1] , rank , node->right->split->kd_coord[0], node->right->split->kd_coord[1]);
      print_tree_P2P(node->right);
    }
  }
}


void debug_tree_P2P(const struct kd_node_t *node)
/*
 * Function to visualize the tree generated by node, mainly for debugging/checking
 */
{
  if(node==NULL) {return;}
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(node->kid==NULL && node->left==NULL && node->right==NULL) {
    printf("{%d}[%d] (%.3f,%.3f) -> NULL \n",
	   node->level, rank, node->split->kd_coord[0], node->split->kd_coord[1] );
    return;
  }
  if(node->kid) {
    if(node->split && node->kid->split) {
      
      printf("{%d}[%d] (%.3f,%.3f) -> (%.3f,%.3f) \n",
	     node->level, rank, node->split->kd_coord[0], node->split->kd_coord[1],
	     node->kid->split->kd_coord[0], node->kid->split->kd_coord[1]);
    }
    debug_tree_P2P(node->kid);
  }
  else {
    if (node->left) {
      printf("{%d}[%d] (%.3f,%.3f) -> (%.3f,%.3f)  \n",
	     node->level, rank, node->split->kd_coord[0], node->split->kd_coord[1],
	     node->left->split->kd_coord[0], node->left->split->kd_coord[1]);
      debug_tree_P2P(node->left);
    }
    
    if (node->right) {
      printf("{%d}[%d] (%.3f,%.3f) -> (%.3f,%.3f) \n",
	     node->level, rank, node->split->kd_coord[0], node->split->kd_coord[1],
	     node->right->split->kd_coord[0], node->right->split->kd_coord[1]);
      debug_tree_P2P(node->right);
    }
  }
}


void clear_P2P(struct kd_node_t *root)
/*
 * This function frees the nodes of the tree
 */
{
  if(root==NULL)
    return;
  if(root->kid) {
    // Note here!!
    // If we are in the P2P parallel region, we have mallocated the node->split
    if(root->split)
      free(root->split);
    clear_shared(root->kid);
  }
  else {
    clear_shared(root->left);
    clear_shared(root->right);
  }
  free(root);
}




void load_from_file_distributed_P2P( int argc,  char ** argv, int *N_ptr, int *N_glo_ptr,  const int ndim, struct kd_point_t ** points_ptr, const int N_p, const int rank)
{
  int rem = 0;
  if(argc != 2) {
    if(rank==0)
      printf( "usage: %s filename\n", argv[0] );
    exit (EXIT_FAILURE);
  }
  MPI_Datatype MPI_kd_point, tmp_type;
  int array_of_blocklengths[] = {ndim};
  MPI_Aint array_of_displacements[] = {0};
  MPI_Datatype array_of_types[] = {MPI_float_t};
  MPI_Aint lb, extent;
  MPI_Type_create_struct(1, array_of_blocklengths, array_of_displacements, array_of_types, &tmp_type);
  MPI_Type_get_extent(tmp_type, &lb, &extent);
  MPI_Type_create_resized(tmp_type, lb, extent, &MPI_kd_point);
  MPI_Type_commit(&MPI_kd_point);
  MPI_Status status;
  int tag = 0;
  FILE *fptr;
  fptr = NULL;
  if(rank == 0) {
    fptr = fopen(argv[1],"r");
    if(fptr == NULL) {
      fprintf(stderr, "%s", "Error occurred while opening file!\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      exit(1);
    }
    printf("MPI kd-tree\n");
    printf("Loadind data from file: %s\n", argv[1]);
    fscanf(fptr,"%d", N_glo_ptr);
    rem = *N_glo_ptr % N_p;
    *N_ptr = *N_glo_ptr/N_p;
  }
  MPI_Bcast(N_ptr, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0) {
    struct kd_point_t one_point;
    for (int n=1; n < N_p; ++n) {
      *points_ptr= (struct kd_point_t *) malloc((*N_ptr)*sizeof(struct kd_point_t));
      for (int i = 0; i < (*N_ptr) ; ++i) {
	for(int d = 0; d < ndim; d++) {
#if defined(DOUBLE_PRECISION)
	  fscanf(fptr,"%lf", &one_point.kd_coord[d]);
#else
	  fscanf(fptr,"%f", &one_point.kd_coord[d]);
#endif
	}
	(*points_ptr)[i] = one_point;
      }
      MPI_Ssend((*points_ptr), *N_ptr, MPI_kd_point, n, tag, MPI_COMM_WORLD);
      free((*points_ptr) );
    }
  }
  else {
    (*points_ptr) = (struct kd_point_t *) calloc((*N_ptr),sizeof(struct kd_point_t));
    MPI_Recv( (*points_ptr), (*N_ptr), MPI_kd_point, 0, tag, MPI_COMM_WORLD, &status);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0) {
    struct kd_point_t one_point;
    (*N_ptr)+= rem ;
    (*points_ptr) = (struct kd_point_t *) calloc((*N_ptr),sizeof(struct kd_point_t));
    for (int i = 0; i < (*N_ptr) ; ++i) {
      for(int d = 0; d< ndim; d++) {
#if defined(DOUBLE_PRECISION)
	fscanf(fptr,"%lf", &one_point.kd_coord[d]);
#else
	fscanf(fptr,"%f", &one_point.kd_coord[d]);
#endif
      }
      (*points_ptr)[i] = one_point;
    }
    fclose(fptr);
  }
  
  MPI_Type_free(&MPI_kd_point);
  MPI_Type_free(&tmp_type);
}



void generate_random_points_distributed_P2P( int argc,  char ** argv, int *N_ptr, int *N_glo_ptr, const int ndim, struct kd_point_t ** points_ptr, const int N_p, const  int rank)
{
  srand48(SEED*(rank+100)) ;
  int rem = 0;
  if(argc != 2) {
    if(rank==0)
      printf( "usage: %s N\n", argv[0] );
    exit (EXIT_FAILURE);
  }
  
  if(argc==2) ( *N_glo_ptr)=atoi(argv[1]);
  
  (*N_ptr) = (*N_glo_ptr) / N_p;
  rem = (*N_glo_ptr) % N_p;
  
  if(rank==0) {
    printf("MPI kd-tree\n");
    printf("Generating %d random double in [0,1]\n", *N_glo_ptr);
  }
  
  if(rank==0) (*N_ptr) += rem;
  
  *points_ptr = (struct kd_point_t*) malloc((*N_ptr)* sizeof(struct kd_point_t));
  for (int i = 0; i < (*N_ptr); i++) {
    for (int j = 0; j < ndim; j++)
      (*points_ptr)[i].kd_coord[j]=drand48();
  }  
}

