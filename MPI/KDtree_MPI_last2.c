#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

// ==================================================================
//                      INPUT parameters
// ==================================================================
//#define DATA "../test_data.csv"

#ifndef NPS
#define NPS 1000
#endif 

#define NDIM 2

//#define FILE_INP "../test_data.csv"

#define ENDLINE "\n"

//#define DOUBLE_PRECISION

//#define DEBUG

// ==================================================================
//                      DATA types and format
// ==================================================================


#ifndef DOUBLE_PRECISION
#define data_t float
#define MPI_DATA_T MPI_FLOAT
#define FMT "%f"
#define PFMT "%6.3f "
#else
#define data_t double
#define MPI_DATA_T MPI_DOUBLE
#define FMT "%lf"
#define PFMT "%6.3lf "
#endif


// ==================================================================
//                       NODE definition
// ==================================================================
typedef struct node KDnode_t;
struct node
{
    data_t coords_point[NDIM];
    int axis, left, right;
    int depth;
};

// declaring MPI datatype as global variable
MPI_Datatype MPI_KDnode_T;

// ==================================================================
//                        FUNCTIONS 
// ==================================================================


/* * * * * * * * * * * * * * * * * * * * * * * * * * *
 Function: Read_File
_______________________________________________________
* Action: reads input file and stores points into an
* array. 
* @param : mpi_rank -> rank of MPI process
* @return : dataset as array of nodes.
* * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifdef FILE_INP
KDnode_t *Read_File(int mpi_rank)
{
    // allocate nodes
    if (mpi_rank == 0)
    {
        printf("N points:%d \n", NPS);
        KDnode_t *tree = (KDnode_t *)malloc(NPS * sizeof(KDnode_t));
  

        // open input file
        char filename[50] = FILE_INP;
        
        FILE *fp = fopen(filename, "r");
        if (fp == NULL)
        {
            perror("Unable to open input file! Exiting...\n");
            exit(1);
        }
        else
        {
            printf("input file opened\n");
        }

        // set format for reading
        char data_fmtsep[4] = FMT;
        char data_fmtlast[6] = FMT;
        char sepstr[2] = ",";
        char endline[2] = ENDLINE;
       
        strncat(data_fmtsep, sepstr, 1);
        strncat(data_fmtlast, endline, 1);
        //printf("%s\n", sepstr);
        //strncat(FMTSEP, sepstr, 1);
 
        
        //printf("I am here\n");

        int r;
        for (int np = 0; np < NPS; ++np)
        {
            for (int nc = 0; nc < NDIM - 1; ++nc)
            {
                r = fscanf(fp, data_fmtsep, (tree + np)->coords_point + nc);
                if (r == EOF)
                {
                    perror("Unexpected end of file while reading!!");
                    exit(2);
                }
            }
            r = fscanf(fp, data_fmtlast, (tree + np)->coords_point + NDIM - 1);
            if (r == EOF)
            {
                perror("Unexpected end of file while reading!!");
                exit(3);
            }
        }

        
        fclose(fp);
        return tree;
    }
    else
        return NULL;
}
#endif

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* Function: Random_Gen
*---------------------------------------------------------------
* @param: mpi_rank, rank of MPI process
*
* returns: pointer to array of points randomly distributed
           for root MPI process (rank 0)

           NULL for other MPI process
 
________________________________________________________________
* DESCRIPTION
  Generates array of points randomly distributed, to use as input
  dataset
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
KDnode_t *Random_Gen(int mpi_rank)
{
  if (mpi_rank==0)
  {
    KDnode_t *tree = (KDnode_t*)malloc(NPS*sizeof(KDnode_t));
  
    srand(1000);

    for (int i = 0; i < NPS; ++i)
    {
           for (int j = 0; j < NDIM; ++j)
           {
                (tree+i)->coords_point[j] = rand()/((double)RAND_MAX)*100;
           }
     
    }

    return tree;

  }


  else {
    return NULL;

  }
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* Function: Check_Tree
*------------------------------------------------------------
* @param: *tree, pointer to the first node in the tree
* @param: mpi_rank, rank of MPI process
*
_____________________________________________________________
  DESCRIPTION
* Prints the tree, as a table, to standard output.
* Table columns are:
* -depth  : depth of KDtree
* -IDX    : index of KDnode in the tree array
* -coords : array of coords of point in KDnode
* -axis   : splitting dimension at next depth
* -IDX_cn : indices of KDnode children in the tree array,
*           when these indices are -1 it means no children
*           and if both the indices are -1, that KDnode is a LEAF
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void Check_Tree(KDnode_t *tree, int mpi_rank)
{
    
    if (mpi_rank == 0)
    {
        //char coordi_format[8] = PFMT; //printing format of not-last coord in data point
        //char coordl_format[8] = PFMT; //printing format of last coord in data point
        //char sepcoord[2] = ",";
        //char endcoord[2] = ")";
        //strncat(coordi_format, sepcoord, 1);
        //strncat(coordl_format, endcoord, 1);
      
        printf("\n");
        printf("IDX   |  KD-NODE coords   |  axis | ID_c_L, ID_c_R | \n");
        for (int np = 0; np < NPS; ++np)
        {
            printf(" %3d | ( ", np);
            for (size_t nc = 0; nc < NDIM; ++nc)
            {
                
                if (nc < NDIM - 1) {
                  printf(PFMT, (tree + np)->coords_point[nc]);
                  printf(",");
                }
                else {
                  printf(PFMT, (tree + np)->coords_point[nc]);
                  printf(")");

                }
            }
            printf(" |     %2d |  %3d  ,  %3d |\n",
                   (tree + np)->axis, (tree + np)->left, (tree + np)->right);
        }
        printf("\n");
    }
}


#ifdef DEBUG
void Check_Tree_mod(KDnode_t *tree)
{
    if (NDIM == 2)
       printf("depth | idx | coord1  | coord2  | axis | id_c_left | id_c_right | \n");
    else if (NDIM == 3)
       printf("depth | idx | coord1 | coord2 | coord3 | axis | id_c_left | id_c_right | \n");
    for (int np = 0; np < NPS; ++np)
    {
        
        printf("%3d   |", (tree+np)->depth);
        printf("%3d  | ", np);
        for (size_t nc = 0; nc < NDIM; ++nc)
        {
            if (nc < (NDIM - 1)) {
               printf(PFMT, (tree + np)->coords_point[nc]);
               printf(" | ");
            }
            else
               printf(PFMT, (tree + np)->coords_point[nc]);
        }
        printf(" | %2d   |    %5d    |    %5d     |\n",
               (tree + np)->axis, (tree + np)->left, (tree + np)->right);
    }
}
#endif

inline void Swap(KDnode_t *a, KDnode_t *b)
{
    data_t tmp[NDIM];
    memcpy(tmp, a->coords_point, sizeof(tmp));
    memcpy(a->coords_point, b->coords_point, sizeof(tmp));
    memcpy(b->coords_point, tmp, sizeof(tmp));
}

int Find_Median(KDnode_t *data, int start, int end, int axis)
{
    if (end <= start)
        return -1;
    if (end == start + 1)
        return start;

    int p, store;
    int m_id = start + (end - start) / 2;
    double pivot;

    while (1)
    {
        // take median as pivot
        pivot = (data + m_id)->coords_point[axis];

        // swap median with end
        Swap((data + m_id), (data + end - 1));

        for (p = store = start; p < end; ++p)
        {
            if ((data + p)->coords_point[axis] < pivot)
            {
                if (p != store)
                    Swap((data + p), (data + store));
                ++store;
            }
        }

        Swap((data + store), (data + end - 1));

        if ((data + store)->coords_point[axis] == (data + m_id)->coords_point[axis])
            return m_id;

        if (store > m_id)
            end = store;
        else
            start = store;
    }
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* Function: Split_Comm
* ----------------------------------------------------------------
 @param: comm, original MPI Communicator
 @param: new_comm, new MPI Communicator
__________________________________________________________________
* DESCRIPTION
* It splits the global MPI communicator in two sub-communicators
* dividing odd and even ranks
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void Split_Comm(MPI_Comm comm, MPI_Comm *new_comm)
{
    int rank, color, key;

    MPI_Comm_rank(comm, &rank);
    color = rank % 2;
    key = rank / 2;

    MPI_Comm_split(comm, color, key, new_comm);
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* Recursive Function: Grow_Tree_Serial
*------------------------------------------------------------------
  @param: *tree, pointer to the KDtree node
  @param: start, index of first point in the branch
  @param: end, index of last point in the branch
  @param: offset, offset needed to store the node in the KDtree
  @param: axis, splitting axis (it can be 0 or 1, for NDIMS=2)

  @returns: m_id+offset, index of new node in the KDtree
            sum of the median point in the branch and offset
___________________________________________________________________
* DESCRIPTION
* It grows the tree serially, up to the point where
* it receives single points (end==start)
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int Grow_Tree_Serial(KDnode_t *tree, int start, int end, int offset, int depth, int axis)
{
    int next_depth = depth + 1;
    // when len of input data is 0 return -1
    if (!(end-start)) return -1; //this case identify the leaves

    // else do the recursive procedure
    int m_id;
    if ((m_id = Find_Median(tree, start, end, axis)) >= 0 ) 
    {
       #ifdef DEBUG
        printf(" at depth %d \n", depth);
       #endif
        (tree+m_id)->depth = next_depth - 1; 
        (tree+m_id)->axis = axis;     //actual splitting axis
        int new_axis = (axis+1) % NDIM; //new is chosen via ROUND-ROBIN

        (tree+m_id)->left  = Grow_Tree_Serial(tree, start, m_id, offset, next_depth, new_axis);
        (tree+m_id)->right = Grow_Tree_Serial(tree, m_id+1, end, offset, next_depth, new_axis);
    }
    return m_id+offset;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* Recursive Function: Grow_Tree_Par
*---------------------------------------------------------------
  @param: *tree, pointer to the KDtree node
  @param: start, index of first point in the branch
  @param: end , index of last point
  @param: offset, offset needed to store the node in the KDtree
  @param: axis, splitting axis (it can be 0 or 1, for NDIMS=2)
  @param: comm, MPI Communicator

  @return: m_id+offset, index of new node in the KDtree
            sum of the median point in the branch and offset
________________________________________________________________
* DESCRIPTION
* It grows the tree in parallel using MPI_Comm_Split,
* to divide the MPI processes in pair of comms
* and assigns the left branches of the tree to
* even rank processes and the right branches to
* odd ranks.
* It calls the function to build the KDtree serially
* when we finishes the MPI processes, meaning 
* when the size of the last communicator is 1
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int Grow_Tree_Par(KDnode_t *tree, int start, int end, int offset, int depth, int axis, MPI_Comm comm)
{
    // get size and ranks in communicator
    int comm_size, comm_rank;
    MPI_Comm_size(comm, &comm_size);
    MPI_Comm_rank(comm, &comm_rank);

    // single process in communicator -> it means that we have finished the MPI processes available
    if (comm_size <= 1)
    {
        return Grow_Tree_Serial(tree, start, end, offset, depth, axis);
    }

    // parallel case
    int m_id = 0, right_count = 0;
    int nleft = -1, nright = -1;
    int next_depth = depth + 1;
    MPI_Status status;

    if (comm_rank == 0)
    {
        m_id = Find_Median(tree, start, end, axis);
       #ifdef DEBUG
        printf(" at depth %d \n", depth);
       #endif
        (tree+m_id)->axis = axis;
        (tree+m_id)->depth = next_depth - 1;
    }
    
    MPI_Bcast(&m_id, 1, MPI_INT, 0, comm);
    right_count = (end - start) - (m_id + 1);

    if (comm_rank == 0)
    {
        // find median and amount of data to send
        MPI_Send((tree+m_id+1), right_count, MPI_KDnode_T, 1, 10*comm_size, comm);
    }
    else if (comm_rank == 1)
    {
        //we allocate the memory for the right branch
        tree = (KDnode_t *)malloc(right_count*sizeof(KDnode_t));
        //we receive the right branch from rank 0
        MPI_Recv(tree, right_count, MPI_KDnode_T, 0, 10*comm_size, comm, &status);
    }
    
    // update axis
    axis = (axis + 1) % NDIM;

    // Splitting the communicator 
    MPI_Comm new_comm;
    Split_Comm(comm, &new_comm);

    // call next step with new communicators
    // the even rank processes take care of the left branch
    if (comm_rank % 2 == 0)
    {
        nleft = Grow_Tree_Par(tree, 0, m_id, offset, next_depth, axis, new_comm);
    }
    // the odd rank processes take care of the right branch
    else
    {
        nright =  Grow_Tree_Par(tree, 0, right_count,  offset+m_id+1, next_depth, axis, new_comm);
    }

    //merging sub-trees in rank 0
    if (comm_rank == 0)
    {   
        //first, we receive from rank 1 the index of right node child in the right branch
        MPI_Recv(&nright, 1, MPI_INT, 1, 11*comm_size, comm, &status);
        (tree+m_id) -> left = nleft;
        (tree+m_id) -> right = nright;
        //after, we receive the data points using the derived MPI datatype
        MPI_Recv(tree+m_id+1, right_count, MPI_KDnode_T, 1, 14*comm_size, comm, &status);    
    }
    else if (comm_rank == 1)
    {
        MPI_Send(&nright, 1, MPI_INT, 0, 11*comm_size, comm);
        // send reordered points
        MPI_Send(tree, right_count, MPI_KDnode_T, 0, 14*comm_size, comm);
        free(tree);
    }

    return m_id + offset;
}

// ==================================================================
//                          MAIN PROGRAM
// ==================================================================
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    // define custom MPI datatype for passing kdnode
    int len_block[2]         = {NDIM, 4};
    MPI_Aint disp[2]         = {0, 2*sizeof(data_t)}; //informs MPI processes about addres of data in the struct
    //ALTERNATIVELY we can use the function MPI_Get_address to get the displacement
    //example: 
    //MPI_Aint disp[2];
    //MPI_Aint base_address;
    //struct KDnode_t dummy_node;
    //MPI_Get_address(&dummy_node, &base_address);
    //MPI_Get_address(&dummy_node.coords_point, &disp[0]);
    //MPI_Get_address(&dummy_node.axis, &disp[1]);
    //disp[0] = MPI_Aint_diff(disp[0], base_address);
    //disp[1] = MPI_Aint_diff(disp[1], base_address);
    MPI_Datatype intypes[2] = {MPI_DATA_T, MPI_INT}; //types of MPI datatype in the struct
    MPI_Type_create_struct(2, len_block, disp, intypes, &MPI_KDnode_T);
    MPI_Type_commit(&MPI_KDnode_T);

    // get size and ranks
    int mpi_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    // read file

 #ifdef FILE_INP
    KDnode_t *tree = Read_File(mpi_rank);
 #else
    KDnode_t *tree = Random_Gen(mpi_rank);
 #endif

    if(mpi_rank==0) {
      printf("N data points: %d \n", NPS);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    double time = MPI_Wtime();
    

    // grow the tree 
    int root = Grow_Tree_Par(tree, 0, NPS, 0, 0, 0, MPI_COMM_WORLD);
    time = MPI_Wtime() - time;
    if(mpi_rank==0) {
      printf("time to build the tree: %lf \n", time);
    }

  
 #ifdef DEBUG
    // print tree for debug
    Check_Tree(tree, mpi_rank);
 #endif

    // average runtime between different MPI processes
    double avg_time;
    MPI_Reduce(&time, &avg_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (mpi_rank == 0)
    {
        // print information about the tree
        printf("Tree grown in %lfs\n", avg_time/mpi_size);
        printf("Tree root is at node %d\n\n", root);
    }

    MPI_Finalize();
    return 0;
}
