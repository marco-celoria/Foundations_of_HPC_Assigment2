#if !defined(DOUBLE_PRECISION)
#define MPI_float_t MPI_FLOAT
#else
#define MPI_float_t MPI_DOUBLE
#endif
#define MIN_NUM_OF_POINTS 4
#define MAX_MPI_LEVEL 5
#define SEED 1234567

void build_kdtree_distributed_A2A(struct kd_node_t * node, struct kd_point_t *points,
                                  int N, int ndim, int axis, int level, MPI_Comm comm);
/*
 * Function to build the parallel MPI kd-tree using the AllToAll simple routine
 */

bool IsPowerOfTwo(int x);
/*
 * Check if x is a power of 2
 */

void clear_A2A(struct kd_node_t *root);
/*
 * This function frees the nodes of the tree
 * NOTE! THIS IS DIFFERENT FROM THE P2P VERSION
 */

void print_tree_A2A(const struct kd_node_t *node);
/*
 * Function to visualize the tree generated by node, mainly for debugging/checking
 */

void debug_tree_A2A(const struct kd_node_t *node);


void show_array_A2A(const  kd_point *data,
                    const int start,
                    const int end,
                    const int level,
                    const char c);
/*
 * Function to print the points in data from start to end,
 * according to the level in the tree
 */


void load_from_file_distributed_A2A    ( int argc,  char ** argv, int *N_ptr, int *N_glo_ptr, const int ndim, struct kd_point_t ** points_ptr, const int N_p, const int rank);


void generate_random_points_distributed_A2A ( int argc,  char ** argv, int *N_ptr, int *N_glo_ptr, const int ndim, struct kd_point_t ** points_ptr, const int N_p, const  int rank);



