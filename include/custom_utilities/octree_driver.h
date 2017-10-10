#pragma once

//system files
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <limits.h>

//C++ system files
#include <vector>
#include <fstream>
#include <bitset>
#include <string>
#include <iomanip>
#include <cstdio>

//my files and other files
#include "geometric_operations.h"
#include "octree_binary.h"
#include "octree_binary_cell.h"
#include <mpi.h>

using namespace Kratos;

#define KRATOS_INDEPENDENT

//relation between lineal position of nodes and position of sons. Lineal positions are numbered following GiD criteria for hexas (right hand rule).
const int sons_positions[ 8 ] = { 0 , 1 , 3 , 2 , 4 , 5 , 7 , 6 };
const int OCTREE_MAX_LEVELS = 30;//max number of levels for octree. This is to avoid infinite recursions.
const int ROOT_LEVEL = OCTREE_MAX_LEVELS - 1;
typedef std::size_t key_type;//This is the binary code used in the octree

/**
 * Class Node
 *
 * This class contains the information of one node, it contains the coordinates, the key,
 * rank in which will be processed and the amount of neighbours.
 *
 * @param coords Is the X, Y and Z coordinates.
 * @param keys Is the binary coding corresponding to the coordinates. 
 * @param rank Index of the rank to process the point.
 * @param n_neighbours The number of neighbours found for the point. 
 */
class Node{
  double coords[ 3 ];
  key_type keys[ 3 ];
  int rank;
  int n_neighbours;
  
public:

  /**
   *Contructor for the node class.
   *
   *The constructor requires de coordinates, keys and the rank.
   *
   *@param[in] x_coord Is the X axis coordinate for the point.
   *@param[in] y_coord Is the Y axis coordinate for the point.
   *@param[in] z_coord Is the Z axis coordinate for the point.
   *@param[in] my_rank Is the rank where the point will be processed.
   *@param[in] x_key Is the key associated to the X coordinate
   *@param[in] y_key Is the key associated to the Y coordinate
   *@param[in] z_key Is the key associated to the Z coordinate
   */
  Node( double x_coord , double y_coord , double z_coord , int my_rank , key_type x_key , key_type y_key , key_type z_key ){
    coords[ 0 ] = x_coord;
    coords[ 1 ] = y_coord;
    coords[ 2 ] = z_coord;
    rank = my_rank;
    keys[ 0 ] = x_key;
    keys[ 1 ] = y_key;
    keys[ 2 ] = z_key;
    n_neighbours = 0;
  }

  /**
   *Method to get the coordinate of a point 
   * 
   *Obtains the coordinte of position i_pos from the node
   *
   *@param[in] i_pos Indicates the coordinate to obtain: 0-X, 1-Y and 2-Z coord.
   *@return The coordinate value from position \a i_pos in a double.
   */
  double _GetCoord( int i_pos ){
    return coords[ i_pos ];
  }

  /**
   *Method to get the rank of a node
   *
   *Gets the rank in wich the node is processed
   *
   *@return The rank where the point is processed
   */
  int _GetRank(){
    return rank;
  }

  /**
   *Method to get a key from the node
   *
   *Gets one of the three keys in the point
   *
   *@param[in] i_pos Is the key position to obtain: 0-X,  1-Y and 2-Z position
   *@return The key from position \a i_pos in a \a key_type.
   */
  key_type _GetKey( int i_pos ){
    return keys[ i_pos ];
  }

  /**
   *Obtain the amount of neighbours
   *
   *This method obtains the amount of neighbours for the node
   *
   *@return The amount of neighbours for the node in a int value.
   */
  int _GetNNeighbours(){
    return n_neighbours;
  }
  
  /**
   *Sum neighbour on node
   *
   *This method sum a neighbour in the node.
   */
  void SumNeighbour(  ){
    n_neighbours++;
  }

  /**
   *Node intersects local neighbour node
   *
   *This method calculates if a local neighbour node( neighbour in the same cell) 
   *intersects with the radius of this node.
   *
   *@param[in] radius Is the search radius for the actual node 
   *@param[in] neighbour_node Is the local  neighbour node pointer.
   */
  void Intersects( double radius , Node* neighbour_node ){
    double neighbour_center[ 3 ];
    neighbour_center[ 0 ] = neighbour_node->_GetCoord( 0 );
    neighbour_center[ 1 ] = neighbour_node->_GetCoord( 1 );
    neighbour_center[ 2 ] = neighbour_node->_GetCoord( 2 );
    if(  Distance( coords , neighbour_center ) <= radius  ){
      n_neighbours++;
      //Adding a new neighbour in the local neighbour
      neighbour_node->SumNeighbour();
    }
  }

  /**
   *Node intersects no local neighbour
   *
   *This method calculates if a non local neighbour intersects with the local neighbour 
   *
   *@param[in] radius Is the search radius for the actual node.
   *@param[in] neighbour_coords Is the center of the neighbour non local node.
   */
  void Intersects( double radius , double* neighbour_coords ){
    if(  Distance( coords , neighbour_coords ) <= radius  ){
      n_neighbours++;
    }
  }

};

/**
 *Read points from file
 *
 *This function reads the nodes information from file
 *
 *@param[in] filename Is the name of the file where the points are stored
 *@param[in] refinement_level Is the refinement used to generate the cells for each rank 
 *@param[out] Points This is the vector of pointers to the nodes read from the file.
 */
void ReadPointsInfo( std::string filename , int refinement_level , std::vector<Node*>& Points ){
  FILE *myfile;
  //Opening data file
  const char *name = filename.c_str();
  myfile = fopen( name , "r" );
	if(!myfile) {
		std::cout << "******ERROR: Can not open points file******" << std::endl;
		exit( 0 );
	}
  int n_points , id , rank;
  double coords[ 3 ];
  fscanf(myfile,"%d", &n_points);
  for(  int i_node = 0  ;  i_node < n_points  ;  i_node++  ){
    fscanf( myfile , "%d" , &id ); //myfile >> id;
    fscanf( myfile , "%lf" , &coords[0] ); //myfile >> coords[ 0 ];
    fscanf( myfile , "%lf" , &coords[1] ); //myfile >> coords[ 1 ];
    fscanf( myfile , "%lf" , &coords[2] ); //myfile >> coords[ 2 ];
    
    key_type keys[3];
    
    keys[ 0 ] = static_cast<key_type> ( ( 1 << ROOT_LEVEL ) * coords[ 0 ]);
    keys[ 1 ] = static_cast<key_type> ( ( 1 << ROOT_LEVEL ) * coords[ 1 ]);
    keys[ 2 ] = static_cast<key_type> ( ( 1 << ROOT_LEVEL ) * coords[ 2 ]);

    //For the nodes on the unit cube boundary is substracted 1 in order to can determine
    //the rank where the point will be processed 
    int flag[ 3 ] = { 0 , 0 , 0 };
    for(  int i_dim = 0  ; i_dim < 3  ;  i_dim++  ){
      if(  keys[ i_dim ] >= ( 1<<ROOT_LEVEL )  ){
        keys[ i_dim ] -= 1;
        flag[ i_dim ] = 1;
      }
    }

    //Calculating the rank for the node 
    rank = 0;
    for(  int i_level = 0  ;  i_level < refinement_level  ;  i_level++  ){
      rank *= 8;
      rank += ( ((keys[ 0 ] & (1<<(ROOT_LEVEL-(i_level+1))))/(1<<(ROOT_LEVEL-(i_level+1)))) +
                ((keys[ 1 ] & (1<<(ROOT_LEVEL-(i_level+1))))/(1<<(ROOT_LEVEL-(i_level+1))))*2 +
                ((keys[ 2 ] & (1<<(ROOT_LEVEL-(i_level+1))))/(1<<(ROOT_LEVEL-(i_level+1))))*4 );
    }
    rank = pow( 8 , refinement_level ) - rank -1;
    Node* aux = new Node( coords[ 0 ] , coords[ 1 ] , coords[ 2 ] , rank , keys[ 0 ]+flag[ 0 ] , keys[ 1 ]+flag[ 1 ] , keys[ 2 ]+flag[ 2 ] );
    Points.push_back( aux );
    
  }
  fclose(myfile);
}

/**
 *Structure of the data contained in each cel
 *
 *This structure stores the information that each cell requires to perform the neighbours
 *search
 *
 *@param Points Is the vector of the pointers to the nodes contained in the cell.
 *@param NL_indexes Is the vevtor that contains the indexes of the other rank neighbours. 
 */
struct DataInsideOctreeCell
{
  std::vector<Node*> Points;
  std::vector<int> NL_indexes;

public:

  /**
   *Default constructor
   */
  DataInsideOctreeCell( ){

  }

  /**
   *Default destructor
   */
  ~DataInsideOctreeCell(){

  }

  /**
   *Obtain the number of nodes in the cell
   *
   *This method obtains the number of nodes contained in the cell
   *
   *@return The amount of nodes in the cell using a \a size_t value. 
   */
  size_t _GetNNodes(){
    return Points.size();
  }

  /**
   *Number of rank neighbours 
   *
   *This method obtains the number of rank neighbours in the cell
   *
   *@return The amount of neighbour ranks for the cell.
   */
  size_t _GetNNeighbours(){
    return NL_indexes.size();
  }

  /**
   *Rank index of neighbour
   *
   *This method obtains the index of the neighbour rank located in position \a i_pos 
   *
   *@param[in] i_pos Is the position of the index to be returned.
   *@return The index of the neighbour rank.
   */
  int GetIndex( int i_pos ){
    return NL_indexes[ i_pos ];
  }

  /**
   *Adding a new node to the cell
   *
   *This method add a new node to the vector of nodes contained in the cell
   *
   *@param[in] node Is a node pointer that will be added in the cell.
   */
  void SetNode( Node* node ){
    Points.push_back( node );
  }

  /**
   *Gets a node pointer
   *
   *This method gets the node pointer in position \a i_node
   *
   *@param[in] i_node Is the node index that will be returned
   *@return The node pointer located in position \a i_node
   */
  Node* GetNode( size_t i_node ){
    return Points[ i_node ];
  }

  /**
   *Adding a neighbour rank index
   *
   *This method adds a new rank index in the cell if it has not been added before.
   *
   *@param[in] rank Is the neighbour rank to be added
   */
  void SetNeighbour( int rank ){
    bool flag = false;
    //Checking if rank exist in the list
    for(  size_t i_pos = 0  ;  i_pos < NL_indexes.size()  ;  i_pos++   ){
      if(  NL_indexes[ i_pos ] == rank  ){
        flag = true;
        break;
      }
    }
    //if rank does not exist, it is added and the vector is sorted
    if(  !flag  ){
      NL_indexes.push_back( rank );
      sort(NL_indexes.begin(), NL_indexes.end(), [](const int a, const int b) {return a > b; });
    }
  }

  /**
   *Count intersections in the cell
   *
   *This method counts the number of neighbours for a set of nodes contained in the cell
   *
   *@param[in] radius Is the search radius for the node.
   */
  void CountLocalNeighbours( double radius ){
    if(  Points.size() > 1  ){
      size_t n_nodes = Points.size();
      for(  size_t i_node = 0  ;  i_node < ( n_nodes - 1  )  ;  i_node++  ){
        Node* node = Points[ i_node ];
        for(  size_t j_node = ( i_node + 1 )  ;  j_node < n_nodes  ;  j_node++  ){
          Node* other_node = Points[ j_node ];
          node->Intersects( radius , other_node );
        }
      }
    }
  }

  /**
   *Count intersection with neighbour cell
   *
   *This method counts the intersections ocurred between the local nodes and the neighbour 
   *cell nodes
   *
   *@param[in] radius Is the search radius for the nodes
   *@param[in] data If the data pointer of the neighbour cell
   */
  void CountNeighbourCellNeighbours( double radius , DataInsideOctreeCell* data ){
    size_t local_nodes = Points.size();
    size_t neighbour_nodes = data->_GetNNodes();
    for(  size_t i_node = 0  ;  i_node < local_nodes  ;  i_node++  ){
      Node* local_node = Points[ i_node ];
      for(  size_t j_node = 0  ;  j_node < neighbour_nodes  ;  j_node++  ){
        Node* neighbour_node = data->GetNode( j_node );
        local_node->Intersects( radius , neighbour_node );
      }
    } 
  } 

  /**
   *Count intersections with neighbour rank
   *
   *This method counts the intersection ocurred between the local nodes and the nodes  
   *received from a neighbour rank
   *
   *@param[in] rank Is the neighbour rank to be processed 
   *@param[in] radius Is the search radius for the nodes
   *@param[in] Recv_coords Is the reference to the vector of coordinates received
   */
  void CountNeighbourRankNeighbours( int rank , double radius , std::vector<std::vector<double>>& Recv_coords ){
    size_t n_nodes = Points.size();
    for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
      Node* node = Points[ i_node ];
      //The number of points is the third part of the size of the vector because each node 
      //has 3 coordinates
      size_t n_points = Recv_coords[ rank ].size()/3;
      for(  size_t i_point = 0  ;  i_point < n_points  ;  i_point++  ){
        double coord[ 3 ];
        coord[ 0 ] = Recv_coords[ rank ][ i_point * 3 + 0 ];
        coord[ 1 ] = Recv_coords[ rank ][ i_point * 3 + 1 ];
        coord[ 2 ] = Recv_coords[ rank ][ i_point * 3 + 2 ];
        node->Intersects( radius , coord );
      }
    }
  }

};  

class OctreeConfigure {
public:
  enum {
    CHILDREN_NUMBER = 8,
    DIMENSION = 3,
    MAX_LEVEL = OCTREE_MAX_LEVELS,
    MIN_LEVEL = 2 // this cannot be less than 2!!! Pooyan
  };
  typedef DataInsideOctreeCell data_type;
  typedef double* point_coords;

  static data_type* AllocateData() {
    return new data_type;
  }

  static void CopyData(data_type* source, data_type* destination) {
    destination = source;
  }

  static void DeleteData(data_type* data) {
    delete data;
  }

};


typedef OctreeBinaryCell<OctreeConfigure> OctreeCell;
typedef OctreeBinary<OctreeCell> Octree;
typedef std::vector<OctreeCell*> OctreeCell_vector;

/**
 *Class comunicaitor, used when a MPI program is executed
 *
 *This class contains the information needed to execute a program using a octree and MPI
 *parallel computing
 *
 *@param mpi_rank Is the rank of the process
 *@param mpi_size Is the number of total process executed.
 *@param g_level Is the refinement needed in the octree to perform the neighbours search
 */
class Comunicator{

	int mpi_rank;
	int mpi_size;
  int g_level;

	public:
    /**
     *Comunicator constructor
     *
     *The constructor initializes the parallel region and calculates the rank, size and 
     *the minimum refinement level required in the octree
     *
     *@param[in] argc A int value that indicates the amount of parameters in the program
     *@param[in] argv If the list of parameters received by the program
     */
		Comunicator(int argc,char **argv){
			MPI_Init( &argc , &argv );
			MPI_Comm_rank( MPI_COMM_WORLD , &mpi_rank );
			MPI_Comm_size( MPI_COMM_WORLD , &mpi_size );
      g_level = (int)(log( mpi_size )/log( 8 ));
		}

    /**
     *Comunicator destructor
     *
     *This destructor finalize the parallel region 
     */
		~Comunicator(){
      MPI_Finalize();
    } 

    /**
     *Get the rank of the process
     *
     *This method obtains the rank of the process thas is being executed
     *
     *@return The rank of the process using a int.
     */
		int  _GetRank(){
      return mpi_rank;
    }

    /**
     *Get number of process executed
     *
     *This method obtains the total amount op process used in the parallel region
     *
     *@return The number of mpi process executed using a int.
     */
		int  _GetSize(){
      return mpi_size;
    }

    /**
     *Get refinement level
     *
     *This method obtains the refinement level needed to create enought cells for each rank
     *
     *@return The refinement level in the octree using a int
     */
		int  _GetLevel(){
      return g_level;
    }

    /**
     *Syncronize all the process
     *
     *This method uses a barrier to syncronize all the mpi process
     */
		void _Synchronize(){
      MPI_Barrier( MPI_COMM_WORLD );
    }

    /**
     *Get rank from a gicen key
     *
     *This method calculates the rank of a node using its key and the global refinement 
     *level
     *
     *@param[in] keys Is the key of the point that the rank is being calculated
     *@return The nuber of rank using a int.
     */
    int GetRankFromKey( key_type* keys ){
      int rank = 0;
      for(  int i_level = 0  ;  i_level < g_level  ;  i_level++  ){
        rank *= 8;
        rank += ( ((keys[ 0 ] & (1<<(ROOT_LEVEL-(i_level+1))))/(1<<(ROOT_LEVEL-(i_level+1)))) +
                  ((keys[ 1 ] & (1<<(ROOT_LEVEL-(i_level+1))))/(1<<(ROOT_LEVEL-(i_level+1))))*2 +
                  ((keys[ 2 ] & (1<<(ROOT_LEVEL-(i_level+1))))/(1<<(ROOT_LEVEL-(i_level+1))))*4 );
      }
      rank = pow( 8 , g_level ) - rank -1;
      return rank;
    }

    /**
     *Sending cell
     *
     *In this method the master process sends to the slaves the cell that each one will 
     *process
     *
     *@param rank Is the index to the rank that the cell will be sended
     *@param cell Is the cell pointer to can get the cell values to send to the slave.
     */
    void MasterSendCellToSlave( int rank , OctreeCell* cell ){
      bool intersection = cell->GetIntersection();//Intersects or not at less with one point
      char level = cell->GetLevel();//Refinement level in the cell
      key_type keys[ 3 ];
      cell->GetKey( 0 , keys );
     
      //26 is the pack size of the cell information
      void* buff = malloc(26);//THIS IS THE SIZE OF 3 SIZE_T, 1 CHAR AND 1 BOOL
      int position = 0;
      MPI_Pack( &level , 1 , MPI_BYTE , buff , 26 , &position , MPI_COMM_WORLD );
      for(  int i_key = 0  ;  i_key < 3  ;  i_key++  )
        MPI_Pack( &keys[ i_key ] , 8 , MPI_BYTE , buff , 26 , &position , MPI_COMM_WORLD );
      MPI_Pack( &intersection , 1 , MPI_BYTE , buff , 26 , &position , MPI_COMM_WORLD );
      MPI_Send( buff , position , MPI_PACKED , rank , 0 , MPI_COMM_WORLD );
      free( buff );
    }
 
    /**
     *Receiving cell
     *
     *This method receives a cell from the master process and creates a instance of the  
     *class \a OctreeCell
     *
     *@return The pointer to the cell created with the information received from master
     */
    OctreeCell* SlaveReceiveCellFromMaster(){
      void* buff = malloc(26);
      MPI_Recv( buff , 26 , MPI_BYTE , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
      OctreeCell *copy;
      bool intersection;
      char level;
      key_type keys[ 3 ];
      //moving in the buffer received to obtain the information packed
      char* ptr = (char*)buff;
      level = *((char*)(ptr));
      ptr += 1;
      for(  int i_key = 0  ;  i_key < 3  ;  i_key++  ){
        keys[ i_key ] = *((key_type*)(ptr));
        ptr += 8;
      }
      intersection = *((bool*)(ptr));
      ptr +=1;
      //creating the new cell with the information received
      copy = new OctreeCell( level , keys[ 0 ] , keys[ 1 ] , keys[ 2 ] , intersection );
      free( buff );
      return copy;
    }

    /**
     *Sending nodes
     *
     *In this method the master process sends the nodes that belongs to a slave process
     *
     *@param[in] rank Is the index to the rank that the nodes will be sent
     *@param[in] Points Is the list of points to be sended to the rank process, the 
     *points are sorted based on the rank of the point, it was performed to reduce the 
     *time when looking for the point to be sent
     */
    void MasterSendNodesToSlave( int rank , std::vector<Node*>& Points ){
      int n_nodes = 0;
      int begin_position = -1;
      int flag = 0;
      //counting the amount of nodes to send
      for(  size_t i_node = 0  ;  i_node < Points.size()  ;  i_node++  ){
        if(  Points[ i_node ]->_GetRank() == rank  ){
          //the flag helps to know if the points counter has started
          if(  flag == 0  ){
            begin_position = (int)i_node;
            flag = 1;
            n_nodes++;
          }else{
            n_nodes++;
          }
        }
      }
      //sending the number of points to be sent
      MPI_Send( &n_nodes , 4 , MPI_BYTE , rank , 0 , MPI_COMM_WORLD );
      if(  n_nodes > 0  ){
        //if at less one node will be sent, then is created a buffer and the information
        //is packed
        void *buff = malloc(24*n_nodes);
        int position = 0;
        for(  int i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
          double coord[ 3 ];
          for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  )
            coord[ i_dim ] = Points[ begin_position + i_node ]->_GetCoord( i_dim );

          for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  )
            MPI_Pack( &coord[ i_dim ] , 8 , MPI_BYTE , buff , 52*n_nodes , &position , MPI_COMM_WORLD );
        }
        MPI_Send( buff , position , MPI_PACKED , rank , 0 , MPI_COMM_WORLD );
        free( buff );
      }
    }

    /**
     *Receiving nodes
     *
     *This method receives a set of points from the master process
     *
     *@param Points[out] Is the list of points received from master process
     */
    void SlaveReceiveNodesFromMaster( std::vector<Node*>& Points ){
      int n_nodes;
      //receives the information about the number of nodes that will receive from master
      MPI_Recv( &n_nodes , 4 , MPI_BYTE , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
      if(  n_nodes > 0  ){
        void* buff = malloc( 24 * n_nodes );
        MPI_Recv( buff , 24 * n_nodes , MPI_BYTE , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
        char *ptr = (char*)buff;
        for(  int i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
          double coord[ 3 ];
          key_type keys[ 3 ];
          int rank;
          for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
            coord[ i_dim ] = *((double*)(ptr));
            ptr += 8;
          }
          keys[ 0 ] = static_cast<key_type> ( ( 1 << ROOT_LEVEL ) * coord[ 0 ]);
          keys[ 1 ] = static_cast<key_type> ( ( 1 << ROOT_LEVEL ) * coord[ 1 ]);
          keys[ 2 ] = static_cast<key_type> ( ( 1 << ROOT_LEVEL ) * coord[ 2 ]);
          Node* aux = new Node( coord[ 0 ] , coord[ 1 ] , coord[ 2 ] , mpi_rank , keys[ 0 ] , keys[ 1 ] , keys[ 2 ] );
          Points.push_back( aux ); 
        }
        free( buff );
      }
    }

    /**
     *Assign comunication information on cells
     *
     *This method check in all the leaves if the neighbours belongs to other process, and  
     * if this happens then adds the neighbour to the information in the data container of
     * the cell.
     *
     *@param[in] local_root Is the pointer to the root of the cell that is being processed
     * in this rank
     */
    void AssignBoundaryInformation( OctreeCell* local_root ){
      //obtaining all leaves
      OctreeCell_vector leaves;
      OctreeCell_vector cells_stack;
      cells_stack.push_back( local_root );
      while( !cells_stack.empty() ){
        OctreeCell* cell = cells_stack.back();
        cells_stack.pop_back();
        if( cell->HasChildren() ){
          for(  size_t i_child = 0  ;  i_child < 8  ;  i_child++  ){
            cells_stack.push_back( cell->pGetChild( i_child ) );
          }
        }else{
          leaves.push_back( cell );
        }
      }
      size_t n_leaves = leaves.size();
      for(  size_t i_leaf = 0  ;  i_leaf < n_leaves  ;  i_leaf++  ){
        OctreeCell* cell = leaves[ i_leaf ];
        if(  cell->GetIntersection()  ){
          DataInsideOctreeCell* data = cell->pGetData();
          for(  size_t i_dir = 0  ;  i_dir < 18  ;  i_dir++  ){
            key_type neighbour_key[ 3 ];
            if(  cell->GetNeighbourKey( i_dir , neighbour_key )  ){
              int rank = GetRankFromKey( neighbour_key );
              if(  rank != mpi_rank  ){
                data->SetNeighbour( rank );
              }
            }
          }
        }else{
          //in this case the cell does not have intersection with any point, so the data
          // in the cell is deleted
          if(  cell->pGetData()  ){
            cell->DeleteData();
          }
        }
      }
    }

    /**
     *Setting local neighbours
     *
     *This method counts the neighbour nodes of the nodes that belong to the local rank
     *
     *@param[in] local_root Is the pointer to the local root cell.
     *@param[in] radius Is the search radius for the cells contained in the rank
     */
    void SetLocalNeighbours( OctreeCell* local_root , double radius ){
      //Obtaining all leavea vector
      OctreeCell_vector leaves;
      OctreeCell_vector cells_stack;
      cells_stack.push_back( local_root );
      while( !cells_stack.empty() ){
        OctreeCell* cell = cells_stack.back();
        cells_stack.pop_back();
        if( cell->HasChildren() ){
          for(  size_t i_child = 0  ;  i_child < 8  ;  i_child++  ){
            cells_stack.push_back( cell->pGetChild( i_child ) );
          }
        }else{
          leaves.push_back( cell );
        }
      }
      size_t n_leaves = leaves.size();
      for(  size_t i_leaf = 0  ;  i_leaf < n_leaves  ;  i_leaf++  ){
        OctreeCell* cell = leaves[ i_leaf ];
        //if the cell intersects with at less one point then the neighbours are searched
        if(  cell->GetIntersection()  ){
          DataInsideOctreeCell* data = cell->pGetData();
          //counting the neighbours located in the same cell 
          data->CountLocalNeighbours( radius );
          //counting neighbours located in the neighbour cells
          size_t index_direction[ 9 ] = {1,3,5,8,9,12,13,16,17};
          for(  size_t i_dir = 0  ;  i_dir < 9  ;  i_dir++  ){
            key_type neighbour_key[ 3 ];
            if(  cell->GetNeighbourKey( index_direction[ i_dir ] , neighbour_key )  ){
              int rank = GetRankFromKey( neighbour_key );
              //in this method are only calculated the neighbours belonging the same rank
              if(  rank == mpi_rank   ){
                //obtaining the neighbour cell
                OctreeCell* neighbour_cell = local_root;
                for(  size_t i_level = 0  ;  i_level < ROOT_LEVEL  ;  i_level++  ){
                    if(  neighbour_cell->IsLeaf()  ) {
                        break;
                    }
                    neighbour_cell = neighbour_cell->pGetChild( neighbour_key[ 0 ] , neighbour_key[ 1 ] , neighbour_key[ 2 ] );
                }
                //if neighbour has points, then are counted the neighbours
                if(  neighbour_cell->GetIntersection( )  ){
                  DataInsideOctreeCell* neighbour_data = neighbour_cell->pGetData();
                  data->CountNeighbourCellNeighbours( radius , neighbour_data );
                }
              }
            }
          }

        }
      }
    }

    /**
     *Obtaining point to send to other rank
     *
     *This method obtains the coordinates that will be sent to other process
     *
     *@param[in] leaves Is the list of leaf cells contained in the local octree
     *@param[out] coords Is the list of coordinates that will be sent to other process
     */
    void GetCoordsToSend( OctreeCell_vector& leaves , std::vector<std::vector<double>>& coords ){
      size_t n_leaves = leaves.size();
      for(  size_t i_leaf = 0  ;  i_leaf < n_leaves  ;  i_leaf++  ){
        OctreeCell* cell = leaves[ i_leaf ];
        //If cell has intersection with some points then is possible that has point to sen
        if(  cell->GetIntersection()  ){
          DataInsideOctreeCell* data = cell->pGetData();
          size_t n_neighbours = data->_GetNNeighbours();
          //if the cell has neighbour cell that belongs to other rank then the points are 
          //sent to the other rank
          if(  n_neighbours  ){
            for(  size_t i_neighbour = 0  ;  i_neighbour < n_neighbours  ;  i_neighbour++  ){
              int index = data->GetIndex( (int)i_neighbour );
              size_t n_nodes = data->_GetNNodes();
              for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
                coords[ index ].push_back( data->GetNode( i_node )->_GetCoord( 0 ) );
                coords[ index ].push_back( data->GetNode( i_node )->_GetCoord( 1 ) );
                coords[ index ].push_back( data->GetNode( i_node )->_GetCoord( 2 ) );
              }
            }
          }
        }
      }
    }

    /**
     *Sending and receiving nodes
     *
     *This method sends and receives the nodes information from all its neighbours
     *
     *@param[in] Send Is the list of nodes that will be sent to each process.
     *@param[out] Recv Is the list of pointers to the nodes that will be receivwd from each neighbour process
     */
    void SendAndReceiveCoordinates(  std::vector<std::vector<double>>& Send , std::vector<std::vector<double>>& Recv  ){

      MPI_Request* Send_req = new MPI_Request[ mpi_size ];
      MPI_Request* Recv_req = new MPI_Request[ mpi_size ];

      MPI_Status* Send_stat = new MPI_Status[ mpi_size ];
      MPI_Status* Recv_stat = new MPI_Status[ mpi_size ];

      int* Send_size = new int[ mpi_size ];
      int* Recv_size = new int[ mpi_size ];

      for(  int i_rank = 0  ;  i_rank < mpi_size  ;   i_rank++  ){
        Send_size[ i_rank ] = (int)Send[ i_rank ].size();
      }

      MPI_Alltoall( Send_size , 1 , MPI_INT , Recv_size , 1 , MPI_INT , MPI_COMM_WORLD );

      std::vector<double*> Send_Buffers( mpi_size , nullptr );
      std::vector<double*> Recv_Buffers( mpi_size , nullptr );
  
      for(  int i_rank = 0  ;  i_rank < mpi_size  ;  i_rank++  ){
        Send_Buffers[ i_rank ] = new double[ Send_size[ i_rank ] ];
        Recv_Buffers[ i_rank ] = new double[ Recv_size[ i_rank ] ];
        for(  int i_pos = 0  ;  i_pos < Send_size[ i_rank ]  ;  i_pos++  ){
          Send_Buffers[ i_rank ][ i_pos ] = Send[ i_rank ][ i_pos ];
        }
        MPI_Isend( Send_Buffers[ i_rank ] , Send_size[ i_rank ] , MPI_DOUBLE , i_rank , 0 , MPI_COMM_WORLD , &Send_req[ i_rank ] );
        MPI_Irecv( Recv_Buffers[ i_rank ] , Recv_size[ i_rank ] , MPI_DOUBLE , i_rank , 0 , MPI_COMM_WORLD , &Recv_req[ i_rank ] );
      }
      for(  int i_rank = 0  ;  i_rank < mpi_size  ;  i_rank++  ){
        MPI_Wait( &Send_req[ i_rank ] , &Send_stat[ i_rank ] );
        MPI_Wait( &Recv_req[ i_rank ] , &Recv_stat[ i_rank ] );
        for(  int i_pos = 0  ;  i_pos < Recv_size[ i_rank ]  ;  i_pos++  ){
          Recv[ i_rank ].push_back( Recv_Buffers[ i_rank ][ i_pos ] );
        }
        delete[] Send_Buffers[ i_rank ];
        delete[] Recv_Buffers[ i_rank ];
      }     

      delete[] Send_size;
      delete[] Recv_size;

      delete[] Send_req;
      delete[] Recv_req;

      delete[] Send_stat;
      delete[] Recv_stat;

    }

    /**
     *Search on neighbours from other rank
     *
     *This method searches for the neighbours located in otther rank
     *
     *@param[in] leaves Is the list of leaves of the octree
     *@param[in] Recv_coords Is the list of vectors of received coordinates
     *@param[in] radius Is the search radius of the points.
     */
    void SearchInNeighboursReceived( OctreeCell_vector& leaves , std::vector<std::vector<double>>& Recv_coords , double radius ){
      size_t n_leaves = leaves.size();
      for(  size_t i_leaf = 0  ;  i_leaf < n_leaves  ;  i_leaf++  ){
        OctreeCell* cell = leaves[ i_leaf ];
        if(  cell->GetIntersection()  ){
          DataInsideOctreeCell* data = cell->pGetData();
          size_t n_neighbours = data->_GetNNeighbours();
          if(  n_neighbours  ){
            for(  size_t i_neighbour = 0  ;  i_neighbour < n_neighbours  ;  i_neighbour++  ){
              int rank = data->GetIndex( (int)i_neighbour );
              data->CountNeighbourRankNeighbours( rank , radius , Recv_coords );
            }
          }
        }
      }
    }

    /**
     *Setting global neighbours
     *
     *This method counts the neighbour nodes of the nodes that belong to the other mpi process
     *
     *@param[in] local_root Is the local root cell for the mpi process
     *@param[in] radius Is the search radius of the points.
     */
    void SetGlobalNeighbours( OctreeCell* local_root , double radius ){
      std::vector<std::vector<double>> Send_coords ( mpi_size , std::vector<double> (0) );
      std::vector<std::vector<double>> Recv_coords ( mpi_size , std::vector<double> (0) );
      OctreeCell_vector leaves;
      OctreeCell_vector cells_stack;
      cells_stack.push_back( local_root );
      while( !cells_stack.empty() ){
        OctreeCell* cell = cells_stack.back();
        cells_stack.pop_back();
        if( cell->HasChildren() ){
          for(  size_t i_child = 0  ;  i_child < 8  ;  i_child++  ){
            cells_stack.push_back( cell->pGetChild( i_child ) );
          }
        }else{
          leaves.push_back( cell );
        }
      }
      GetCoordsToSend( leaves , Send_coords );
      SendAndReceiveCoordinates(  Send_coords , Recv_coords  );
      SearchInNeighboursReceived( leaves , Recv_coords , radius );
    }

    /**
     *Searching for neighbours 
     *
     *This method calls to the search in the local neighbours and after that in the 
     *neighbours from other rank
     *
     *@param[in] local_root Is the local root cell for the mpi process
     *@param[in] radius Is the search radius of the points.
     */
    void NeighboursSearch( OctreeCell* local_root , double radius ){
      SetLocalNeighbours( local_root , radius );
      SetGlobalNeighbours( local_root , radius );
    }
};

//class GiDOM_Octree. This class is like an interafce to use OctreeBinary function, taking into account the needs and data of the mesher
struct OctreeDriver
{
  Octree* octree_bin;
public:

  OctreeDriver(){
    octree_bin = new Octree();
  }
  ~OctreeDriver(){
    delete octree_bin;
  }

  Octree* GetOctree() const {return octree_bin;}

	bool CheckLocalBalanceMPI();
	bool CheckGlobalBalanceMPI( Comunicator *com );
	int BalanceOctree();

	OctreeCell* GetRoot() {
    return octree_bin->pGetRoot();
  };

	OctreeCell* GetOctreeCell( const double *coord ) const {
    return octree_bin->pGetCellNormalized( coord );
  }

	OctreeCell* GetOctreeCell( Octree::key_type* keys ) const {
    return octree_bin->pGetCell( keys );
  }

  OctreeCell_vector CalcAllLeavesVector(){
    OctreeCell_vector all_leaves;
	  octree_bin->GetAllLeavesVector( all_leaves );
	  return all_leaves;
  }

  void GetLeavesInBoundingBox( const double* coord1 , const double* coord2,
    OctreeCell_vector& leaves) const {
      return octree_bin->GetLeavesInBoundingBoxNormalized( coord1 , coord2 , leaves );
  }

  int SubdivideCell( OctreeCell* cell ) const;

  int RefineOctreeWithSize( const double size ) {
    int fail = octree_bin->RefineWithUniformSizeNormalized( size );
    return fail;
  }

  OctreeCell* GetNeighbour( const OctreeCell* cell , const int idir ) const {
    OctreeCell* ret = octree_bin->pGetNeighbourCell( cell , idir );
    return ret;
  }

  OctreeCell* GetLeftNeighbour( const OctreeCell* cell ) const {
    OctreeCell* ret= octree_bin->pGetLeftCell( cell );
    return ret;
  }

  OctreeCell* GetRightNeighbour( const OctreeCell* cell ) const {
    OctreeCell* ret= octree_bin->pGetRightCell( cell );
    return ret;
  }
  OctreeCell* GetTopNeighbour( const OctreeCell* cell ) const {
    OctreeCell* ret= octree_bin->pGetTopCell( cell );
    return ret;
  }
  OctreeCell* GetBottomNeighbour( const OctreeCell* cell ) const {
    OctreeCell* ret= octree_bin->pGetBottomCell( cell );
    return ret;
  }
  OctreeCell* GetFrontNeighbour( const OctreeCell* cell ) const{
    OctreeCell* ret= octree_bin->pGetFrontCell( cell );
    return ret;
  }
  OctreeCell* GetBackNeighbour( const OctreeCell* cell ) const{
    OctreeCell* ret= octree_bin->pGetBackCell( cell );
    return ret;
  }

	void RefineWithUniformSizeNormalized( const double level ) {
    octree_bin->RefineWithUniformSizeNormalized( level );
  }

	void GetCoordOctreePosition( const OctreeCell* cell , int ipos , double* coord_point ) const{
	  key_type keys[ 3 ];
	  cell->GetKey(  ipos , keys  );
	  for (  int i = 0  ;  i < 3  ;  i++  ){
		  coord_point[ i ] = GetCoordinate( keys[ i ] );
	  }
  }

	double GetCoordinate( Octree::key_type key ) const {
    return GetOctree()->GetCoordinateNormalized( key );
  }

	double CalcSize( const OctreeCell* cell ) const {
    return GetOctree()->CalcSizeNormalized( cell );
  }

	void PrintMeshGiD( Comunicator *com ){
	  int mpi_rank=com->_GetRank();
	  char rank_name[ 50 ];
	  sprintf( rank_name , "%d" , mpi_rank );
	  std::ofstream rOStream;
	  strcat( rank_name , "mesh.msh" );
	  rOStream.open( rank_name );
	  size_t nleaves = 0;
	  OctreeCell_vector leaves;
	  octree_bin->GetAllLeavesVector( leaves );
	  nleaves = (int)( leaves.size() );
	  //std::cout << mpi_rank << " writing " << nleaves << " leaves" << std::endl;
	  rOStream << "MESH \"leaves\" dimension 3 ElemType Hexahedra Nnode 8\n";
	  rOStream << "# color 96 96 96" << std::endl;
	  rOStream << "Coordinates" << std::endl;
	  rOStream << "# node number coordinate_x coordinate_y coordinate_z  " << std::endl;
	  size_t node_index = 1;
	  for(  size_t i = 0  ;  i < nleaves  ;  i++  ) {
		  OctreeCell* cell = leaves[ i ];
		  double min_point[ 3 ];
		  double max_point[ 3 ];
		  GetCoordOctreePosition( cell , 0 , min_point );
		  GetCoordOctreePosition( cell , 6 , max_point );
		  double cell_size = max_point[ 0 ] - min_point[ 0 ];

		  for(  size_t j = 0  ;  j < 2  ;  j++  )
			  for(  size_t k = 0  ;  k < 2  ;  k++  )
				  for(  size_t h = 0  ;  h < 2  ;  h++  ) {
					  rOStream << node_index++ << "  " << min_point[ 0 ] + j * cell_size << "  " << min_point[ 1 ] + k * cell_size << "  " << min_point[ 2 ] + h * cell_size << std::endl;
				  }

	  }
	  rOStream << "end coordinates" << std::endl;
	  rOStream << "Elements" << std::endl;
	  rOStream << "# element node_1 node_2 node_3 material_number" << std::endl;

	  for(  size_t i = 0  ;  i < leaves.size()  ;  i++  ) {
		  if (  leaves[ i ]->pGetData()  )
			  rOStream << i + 1 << "  " << 8 * i + 1 << "  " << 8 * i + 2 << "  " << 8 * i + 4 << "  " << 8 * i + 3 << "  " << 8 * i + 5 << "  " << 8 * i + 6 << "  " << 8 * i + 8 << "  " << 8 * i + 7 << "  " << leaves[ i ]->GetLevel() + 100 << std::endl;
		  else
			  rOStream << i + 1 << "  " << 8 * i + 1 << "  " << 8 * i + 2 << "  " << 8 * i + 4 << "  " << 8 * i + 3 << "  " << 8 * i + 5 << "  " << 8 * i + 6 << "  " << 8 * i + 8 << "  " << 8 * i + 7 << "  " << int(leaves[ i ]->GetLevel()) << std::endl;

	  }
	  rOStream << "end elements" << std::endl;
	  rOStream.close();

  }
};

/**
 *Setting intersections on rank cells
 *
 *This function set true or false, depending on the existence or not the intersections  
 *from the cell with at less one point
 *
 *@param[in] leaves This is the list of cells in the global octree
 *@param[out] leaves In the data inside octree cell is setted true or false depending on the intersection or not intersection with the points
 *@param[in] Points Is the list of nodes to search its neighbours.
 */
void AsignIntersectionOnCells( OctreeCell_vector& leaves , std::vector<Node*>& Points ){
  size_t n_points = Points.size();
  int* flag = new int [ leaves.size() ];
  size_t cont = 0;
  //Setting no intersections on cells before check for intersections on nodes
  for(  size_t i_leaf = 0  ;  i_leaf < leaves.size()  ;  i_leaf++  ){
    flag[ i_leaf ] = 0;
    leaves[ i_leaf ]->SetIntersection( false );
  }
  for( size_t i_node = 0  ;  i_node < n_points  ;  i_node++ ){
    if(  cont == leaves.size()  ){
      break;
    }
    int rank = Points[ i_node ]->_GetRank();
    if( flag[rank] == 0 ){
      leaves[ rank ]->SetIntersection( true );
      flag[ rank ] = 1;
      cont++;
    }
  }
  delete[] flag;
}

/**
 *Cleaning list of leaves and list of points
 *
 *This method cleans the list of leaves and the list of nodes from the master process in 
 *order to can search in the local information.
 *
 *@param[out] levaes Is the list of leaves that will be erased
 *@param[out] Points Is the list of points, at the end this list will only contain the set of nodes belonging to the local rank
 */
void CleanLeavesAndPoints( OctreeCell_vector& leaves , std::vector<Node*>& Points ){
  leaves.clear();
  size_t total_nodes = Points.size();
  size_t n_nodes = 1;
  int actual_rank = Points[ 0 ]->_GetRank();
  int next_rank;
  for(  size_t i_node = 1  ;  i_node < total_nodes  ;  i_node++  ){
    next_rank = Points[ i_node ]->_GetRank();
    if(  next_rank != actual_rank  ){
      break;
    }
    actual_rank = next_rank;
    n_nodes++;
  }
  size_t i_point = 0;
  while(  i_point < ( total_nodes - n_nodes )  ){
    Points.pop_back();
    i_point++;
  }
}

/**
 *Refine cells intersected by points
 *
 *This function refine the cells intersected by points until reach a refinement level that
 *guarantees that the neighbour nodes are localet in the leave that the node belongs to, or 
 *the neighbours are located on one of the neighbour cells
 *
 *@param[out] local_root Is the pointer to the local root cell of the process
 *@param[in] Points Is the list of points that belongs to the local process
 *@param[in] level is the maximum refinement level that a cell can reach
 */
void RefineCellsIntersected( OctreeCell* local_root , std::vector<Node*>& Points , int level ){

  if( local_root->GetIntersection() ){
    DataInsideOctreeCell* data_ = new DataInsideOctreeCell;
    for(  size_t i_node = 0  ;  i_node < Points.size()  ;  i_node++  ){
      data_->SetNode( Points[ i_node ] );
    }
    DataInsideOctreeCell** data_cell = local_root->pGetDataPointer();
    (*data_cell) = data_;  

    int i_level = 0;
    while( i_level < level ){
      OctreeCell_vector leaves;
      i_level++;
      //Obtaining all leaves vector
      OctreeCell_vector cells_stack;
      cells_stack.push_back( local_root );
      while( !cells_stack.empty() ){
        OctreeCell* cell = cells_stack.back();
        cells_stack.pop_back();
        if( cell->HasChildren() ){
          for(  size_t i_child = 0  ;  i_child < 8  ;  i_child++  ){
            cells_stack.push_back( cell->pGetChild( i_child ) );
          }
        }else{
          leaves.push_back( cell );
        }
      }
		  OctreeCell* leaf;
		  size_t n_leaves = leaves.size();
		  for(  size_t i_leaf = 0  ;  i_leaf < n_leaves  ;  i_leaf++  ){
			  leaf = leaves[ i_leaf ];
        //A cell is subdivided is intersects with at less one point
			  if( leaf->GetIntersection() ){
				  leaf->SubdivideCell(  );
          DataInsideOctreeCell** data_children = new DataInsideOctreeCell*[ 8 ];
          for(  size_t i_child = 0  ;  i_child < 8  ;  i_child++  ){
            data_children[ i_child ] = new DataInsideOctreeCell;
            OctreeCell* child = leaf->pGetChild( i_child );
            DataInsideOctreeCell** data_child = child->pGetDataPointer();
            (*data_child) = data_children[ i_child ]; 
          }

          //Assigning data to children
          DataInsideOctreeCell* data_leaf = leaf->pGetData();
          size_t n_nodes = data_leaf->_GetNNodes();
          for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
            key_type keys[ 3 ];
            Node* node = data_leaf->GetNode( i_node ); 
            keys[ 0 ] = node->_GetKey( 0 );
            keys[ 1 ] = node->_GetKey( 1 );
            keys[ 2 ] = node->_GetKey( 2 ); 
            size_t index = leaf->GetChildIndex(  keys[ 0 ] , keys[ 1 ] , keys[ 2 ]  );
            data_children[ index ]->SetNode( node );
          }
          for(  size_t i_child = 0  ;  i_child < 8  ;  i_child++  ){
            OctreeCell* child = leaf->pGetChild( i_child );
            if(  data_children[ i_child ]->_GetNNodes()  ){
              child->SetIntersection( true );
            }else{
              child->SetIntersection( false );
            }
          }
			  }
		  }
		  leaves.clear();
    }
  }
}

/**
 *Search the neighbours inside radius from a cloud of points
 *
 *  This function search all the neighbours inside a radius from a set of points.
 *
 *  @param filename Is the name of the file containing the cloud of points.
 *  @param radius Is the search radius
 *  @param arg Is the number of parameters recieved in the program
 *  @param argv Is the list of params recieved
 */

void RunMultiplePointMPISearchOctree( std::string filename , double radius , int arg , char* argv[] ){
  
  //Beginning the mpi paralle region
	Comunicator *com=new Comunicator( arg , argv );

  //cell pointer used in each cell
  OctreeCell* local_root = NULL;

  //calculating the maximum refinement level used in the local octree, cell size have to 
  //be the smallest possible but not smaller than the search radius
  int level = com->_GetLevel();
  double min_cell_size = 1.0/pow( 2 , level );
  while( 1 ){
    double next_size = 1.0/pow( 2 , level+1 );
    if( next_size <= radius ){
      break;
    }else{
      level++;
    }
  }
  level = level - com->_GetLevel();


  if( com->_GetRank() == 0 ){//This is the master process

    double t0, t1 , t2; 
    t0 = MPI_Wtime();
    t1 = MPI_Wtime();
    t2 = MPI_Wtime() - t0;
    fprintf(stdout, "[%16lf]-[%16lf]   **MASTER-BEGINNING MULTIPLE POINTS SEARCH\n\n\n", t2 , MPI_Wtime() - t1 );
    
    //Creating octree and refine until have a cell for each mpi_process
    t1 = MPI_Wtime();
    double cell_size = 1.0 / pow( 2 , com->_GetLevel() );
	  OctreeDriver octree_driver;
	  octree_driver.RefineWithUniformSizeNormalized( cell_size );
    t2 = MPI_Wtime() - t0;
		fprintf(stdout, "[%16lf]-[%16lf]   **Refinement             : \n", t2 , MPI_Wtime() - t1 );

    //Obtainig all leaves
    t1 = MPI_Wtime();
	  OctreeCell_vector leaves;
	  leaves = octree_driver.CalcAllLeavesVector( );
    t2 = MPI_Wtime() - t0;
		fprintf(stdout, "[%16lf]-[%16lf]   **Calculating all leaves : \n", t2 , MPI_Wtime() - t1 );

    //Reading points information
    t1 = MPI_Wtime();
    std::vector<Node*> Points;
    ReadPointsInfo( filename , com->_GetLevel() , Points );
    t2 = MPI_Wtime() - t0;
		fprintf(stdout, "[%16lf]-[%16lf]   **Read points information: \n", t2 , MPI_Wtime() - t1 );

    //Determining leaves intersected by points
    t1 = MPI_Wtime();
    AsignIntersectionOnCells( leaves , Points );
    t2 = MPI_Wtime() - t0;
		fprintf(stdout, "[%16lf]-[%16lf]   **Setting intersections  : \n", t2 , MPI_Wtime() - t1 );

    //Sending leaves to slaves
    t1 = MPI_Wtime();
    local_root = leaves[ 0 ];
    for(  size_t i_leaf = 1 ; i_leaf < leaves.size()  ;  i_leaf++  ){
      com->MasterSendCellToSlave( (int)i_leaf , leaves[ i_leaf ] );
    }
    t2 = MPI_Wtime() - t0;
		fprintf(stdout, "[%16lf]-[%16lf]   **Sending leaves         : \n", t2 , MPI_Wtime() - t1 );

    //Send set of points to each slave
    t1 = MPI_Wtime();
    sort( Points.begin( ), Points.end( ), [ ]( Node* lhs, Node* rhs ){ return lhs->_GetRank() < rhs->_GetRank(); });
    for(  size_t i_leaf = 1 ; i_leaf < leaves.size()  ;  i_leaf++  ){
      com->MasterSendNodesToSlave( (int)i_leaf , Points );
    }
    CleanLeavesAndPoints( leaves , Points );
    t2 = MPI_Wtime() - t0;
		fprintf(stdout, "[%16lf]-[%16lf]   **Sending points         : \n", t2 , MPI_Wtime() - t1 );


    //subdividing local root, if intersect some points 
    t1 = MPI_Wtime(); 
    RefineCellsIntersected( local_root , Points , level );
    t2 = MPI_Wtime() - t0;
		fprintf(stdout, "[%16lf]-[%16lf]   **Refining octree        : \n", t2 , MPI_Wtime() - t1 );

    //Assigning boundary information in each cell
    t1 = MPI_Wtime();
    com->AssignBoundaryInformation( local_root );
    t2 = MPI_Wtime() - t0;
		fprintf(stdout, "[%16lf]-[%16lf]   **Setting boundary info  : \n", t2 , MPI_Wtime() - t1 );

    //Searching for neighbours
    t1 = MPI_Wtime();
    com->NeighboursSearch( local_root , radius );
    t2 = MPI_Wtime() - t0;
		fprintf(stdout, "[%16lf]-[%16lf]   **Neighbours search      : \n", t2 , MPI_Wtime() - t1 );  

    fprintf(stdout, "                                        **_______________________:\n");
    t2 = MPI_Wtime() - t0;
		fprintf(stdout, "[%16lf]                      **TOTAL TIME             : \n", t2 );

  }else{ //Slaves
    

    std::vector<Node*> Points;

    //Receiving information from master
    local_root=com->SlaveReceiveCellFromMaster( );
    com->SlaveReceiveNodesFromMaster( Points );

    //Subdividing local root if intersects some points
    RefineCellsIntersected( local_root , Points , level );

    //Adding boundary information in each cell
    com->AssignBoundaryInformation( local_root );

    //Searching for neighbours
    com->NeighboursSearch( local_root , radius );

  }
	MPI_Finalize();
}






