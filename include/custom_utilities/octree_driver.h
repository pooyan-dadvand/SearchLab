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

//my files and other files
#include "octree_binary.h"
#include "octree_binary_cell.h"
#include "mpi.h"

using namespace Kratos;

#define KRATOS_INDEPENDENT

//relation between lineal position of nodes and position of sons. Lineal positions are numbered following GiD criteria for hexas (right hand rule).
const int sons_positions[ 8 ] = { 0 , 1 , 3 , 2 , 4 , 5 , 7 , 6 };
const int OCTREE_MAX_LEVELS = 30;//max number of levels for octree. This is to avoid infinite recursions.
const int ROOT_LEVEL = OCTREE_MAX_LEVELS - 1;
typedef std::size_t key_type;

class Node{
  double coords[ 3 ];
  key_type keys[ 3 ];
  int rank;
public:
  Node( double x_coord , double y_coord , double z_coord , int my_rank , key_type x_key , key_type y_key , key_type z_key ){
    coords[ 0 ] = x_coord;
    coords[ 1 ] = y_coord;
    coords[ 2 ] = z_coord;
    rank = my_rank;
    keys[ 0 ] = x_key;
    keys[ 1 ] = y_key;
    keys[ 2 ] = z_key;
  }
  double _GetCoord( int i_pos ){
    return coords[ i_pos ];
  }
  int _GetRank(){
    return rank;
  }
  key_type _GetKey( int i_pos ){
    return keys[ i_pos ];
  }
};

std::vector<Node> ReadPointsInfo( std::string filename , int refinement_level ){
  std::ifstream myfile;
  //OPENING POINTS DATA FILE
  myfile.open( filename );
  std::vector<Node> Points;
  int n_points , trash , rank;
  double coords[ 3 ];
  myfile >> n_points;
  for(  int i_node = 0  ;  i_node < n_points  ;  i_node++  ){
    myfile >> trash;
    myfile >> coords[ 0 ];
    myfile >> coords[ 1 ];
    myfile >> coords[ 2 ];
    
    key_type keys[3];
    
    keys[ 0 ] = static_cast<key_type> ( ( 1 << ROOT_LEVEL ) * coords[ 0 ]);
    keys[ 1 ] = static_cast<key_type> ( ( 1 << ROOT_LEVEL ) * coords[ 1 ]);
    keys[ 2 ] = static_cast<key_type> ( ( 1 << ROOT_LEVEL ) * coords[ 2 ]);

    int flag[ 3 ] = { 0 , 0 , 0 };
    if(  keys[ 0 ] >= ( 1<<ROOT_LEVEL )  ){
      keys[ 0 ] -= 1;
      flag[ 0 ] = 1;
    }
    if(  keys[ 1 ] >= ( 1<<ROOT_LEVEL )  ){
      keys[ 1 ] -= 1;
      flag[ 1 ] = 1;
    }
    if(  keys[ 2 ] >= ( 1<<ROOT_LEVEL )  ){
      keys[ 2 ] -= 1;
      flag[ 2 ] = 1;
    }

    rank = 0;
    for(  int i_level = 0  ;  i_level < refinement_level  ;  i_level++  ){
      rank *= 8;
      rank += ( ((keys[ 0 ] & (1<<(ROOT_LEVEL-(i_level+1))))/(1<<(ROOT_LEVEL-(i_level+1)))) +
                ((keys[ 1 ] & (1<<(ROOT_LEVEL-(i_level+1))))/(1<<(ROOT_LEVEL-(i_level+1))))*2 +
                ((keys[ 2 ] & (1<<(ROOT_LEVEL-(i_level+1))))/(1<<(ROOT_LEVEL-(i_level+1))))*4 );
    }
    rank = pow( 8 , refinement_level ) - rank -1;
    Node aux( coords[ 0 ] , coords[ 1 ] , coords[ 2 ] , rank , keys[ 0 ]+flag[ 0 ] , keys[ 1 ]+flag[ 1 ] , keys[ 2 ]+flag[ 2 ] );
    Points.push_back( aux );

  }  
  //CLOSING POINTS DATA FILE
  myfile.close();
  return Points;
}

struct DataInsideOctreeCell
{
  std::vector<Node> Points;
public:
  DataInsideOctreeCell(){
    
  }
  ~DataInsideOctreeCell(){

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


class Comunicator{

	int mpi_rank;
	int mpi_size;

	public:
		Comunicator(int argc,char **argv){
			MPI_Init( &argc , &argv );
			MPI_Comm_rank( MPI_COMM_WORLD , &mpi_rank );
			MPI_Comm_size( MPI_COMM_WORLD , &mpi_size );
		}

		~Comunicator(){
      MPI_Finalize();
    }

		int  _GetRank(){
      return mpi_rank;
    }

		int  _GetSize(){
      return mpi_size;
    }

		void _Synchronize(){
      MPI_Barrier( MPI_COMM_WORLD );
    }

    void MasterSendCellToSlave( int rank , OctreeCell* cell ){
      bool intersection = cell->GetIntersection();
      char level = cell->GetLevel();
      key_type keys[ 3 ];
      cell->GetKey( 0 , keys );
     
      void* buff = malloc(26);//THIS IS THE SIZE OF 3 SIZE_T, 1 CHAR AND 1 BOOL
      int position = 0;
      MPI_Pack( &level , 1 , MPI_BYTE , buff , 26 , &position , MPI_COMM_WORLD );
      for(  int i_key = 0  ;  i_key < 3  ;  i_key++  )
        MPI_Pack( &keys[ i_key ] , 8 , MPI_BYTE , buff , 26 , &position , MPI_COMM_WORLD );
      MPI_Pack( &intersection , 1 , MPI_BYTE , buff , 26 , &position , MPI_COMM_WORLD );
      MPI_Send( buff , position , MPI_PACKED , rank , 0 , MPI_COMM_WORLD );
      free( buff );
    }
 
    OctreeCell* SlaveReceiveCellFromMaster(){
      void* buff = malloc(26);
      MPI_Recv( buff , 26 , MPI_BYTE , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
      OctreeCell *copy;
      bool intersection;
      char level;
      key_type keys[ 3 ];
      char* ptr = (char*)buff;
      level = *((char*)(ptr));
      ptr += 1;
      for(  int i_key = 0  ;  i_key < 3  ;  i_key++  ){
        keys[ i_key ] = *((key_type*)(ptr));
        ptr += 8;
      }
      intersection = *((bool*)(ptr));
      ptr +=1;
      copy = new OctreeCell( level , keys[ 0 ] , keys[ 1 ] , keys[ 2 ] , intersection );
      free( buff );
      return copy;
    }

    void MasterSendNodesToSlave( int rank , std::vector<Node> Points ){
      int n_nodes = 0;
      int begin_position = -1;
      int flag = 0;
      for(  size_t i_node = 0  ;  i_node < Points.size()  ;  i_node++  ){
        if(  Points[ i_node ]._GetRank() == rank  ){
          if(  flag == 0  ){
            begin_position = (int)i_node;
            flag = 1;
            n_nodes++;
          }else{
            n_nodes++;
          }
        }
      }
      MPI_Send( &n_nodes , 4 , MPI_BYTE , rank , 0 , MPI_COMM_WORLD );
      void *buff = malloc(52*n_nodes);
      int position = 0;
      for(  int i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
        double coord[ 3 ];
        key_type keys[ 3 ];
        for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
          coord[ i_dim ] = Points[ begin_position + i_node ]._GetCoord( i_dim ); 
          keys[ i_dim ] = Points[ begin_position + i_node ]._GetKey( i_dim ); 
        }
        for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  )
          MPI_Pack( &coord[ i_dim ] , 8 , MPI_BYTE , buff , 52*n_nodes , &position , MPI_COMM_WORLD );

        MPI_Pack( &rank , 4 , MPI_BYTE , buff , 52*n_nodes , &position , MPI_COMM_WORLD );

        for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  )
          MPI_Pack( &keys[ i_dim ] , 8 , MPI_BYTE , buff , 52*n_nodes , &position , MPI_COMM_WORLD );

      }
      MPI_Send( buff , position , MPI_PACKED , rank , 0 , MPI_COMM_WORLD );
      free( buff );
    }

    std::vector<Node> SlaveReceiveNodesFromMaster( ){
      int n_nodes;
      MPI_Recv( &n_nodes , 4 , MPI_BYTE , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
      void* buff = malloc( 52 * n_nodes );
      MPI_Recv( buff , 52 * n_nodes , MPI_BYTE , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
      std::vector<Node> Points;
      char *ptr = (char*)buff;
      for(  int i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
        double coord[ 3 ];
        key_type keys[ 3 ];
        int rank;
        for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
          coord[ i_dim ] = *((double*)(ptr));
          ptr += 8;
        }
        rank = *((int*)(ptr));
        ptr += 4;
        for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
          keys[ i_dim ] = *((size_t*)(ptr));
          ptr += 8;
        }
        Node aux( coord[ 0 ] , coord[ 1 ] , coord[ 2 ] , rank , keys[ 0 ] , keys[ 1 ] , keys[ 2 ] );
        Points.push_back( aux ); 
      }
      free( buff );
      return Points;
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
	  std::cout << mpi_rank << " writing " << nleaves << " leaves" << std::endl;
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
	void AsignRankLeaves( OctreeCell_vector leaves);
	int BalanceOctreeMPI( Comunicator *com );
	void BalanceOctreeIntraprocess();
	void BalanceOctreeInterprocess( Comunicator *com );
	size_t GetRankLeaf( size_t *keys );
	void CopyKeys( int position , size_t *dest , size_t *source );
	bool DoesExistCellWithKey( size_t *test_key , OctreeCell_vector leaves );
	bool CompareKeys( size_t *current_key , size_t *test_key );
	void GetParentKey( size_t *current_key , char current_level );
	int GetChildId( size_t *keys , char level );
	size_t GetNLeavesInsideKeyAndSize( size_t *keys , char level , OctreeCell_vector& all_leaves );
	void GetAllKeysAndLevelsFromLeaves( size_t *keys , char *level , OctreeCell_vector leaves );
	void BalanceLocalCellsBasedOnGhostCells( size_t n_ghost_cells , size_t *ghost_keys , char *ghost_levels , OctreeCell_vector leaves , int my_rank );
	void GetAndSendCellWhichContainsKey( size_t *aux_key , OctreeCell_vector leaves , int process , Comunicator *com );
	bool IsCurrentLeafMinimumThanNeighbourLeaf( size_t* current_key , size_t * neighbour_key );
  key_type GetKeyNormalized( double coord ){
    return octree_bin->CalcKeyNormalized( coord );
  }
};

bool Intersects( double *min_coord , double *max_coord , double *center , double radius  ){
  double dist = 0.0;
  for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
    double edge = center[ i_dim ] - min_coord[ i_dim ];
    if(  edge < 0.0  ){
      if(  edge < -radius  ){
        return false;
      }else{
        dist += ( edge * edge );
      }
    }else{
      edge = center[ i_dim ] - max_coord[ i_dim ];;
      if(  edge > 0.0  ){
        if(  edge > radius  ){
          return false; 
        }else{
          dist += ( edge * edge );
        }
      }
    }
  }
  if(  dist <= ( radius * radius )  ){
    return true;
  }
  return false;
}

void RunPointSearchOctree( std::string filename , double radius , double *center , int level , int arg , char* argv[] ){
  //BEGINING PARALLEL REGION
	Comunicator *com=new Comunicator( arg , argv );

  //THIS IS THE CELL POINTER THAT WILL BE PROCESSED ON EACH MPI PARTITION
  OctreeCell* local_root = NULL;

  if( com->_GetRank() == 0 ){ //THIS IS THE MASTER PROCESS


	  OctreeDriver octree_driver;
    //CELL SIZE CALCULATED TO CREATE CELLS FOR EACH PROCESS
    int refinement_level = log( (double)com->_GetSize() )/log( 8 );
    double cell_size = 1.0 / pow( 2 , refinement_level );
	  octree_driver.RefineWithUniformSizeNormalized( cell_size );

    //READING POINTS POSITION
    std::vector<Node> Points;
    Points = ReadPointsInfo( filename , refinement_level );

    //OBTAINING ALL LEAVES
	  OctreeCell_vector leaves;
	  leaves = octree_driver.CalcAllLeavesVector( );
    
    //CALCULATING INTERSECTIONS
    for (  int i_leaf = 0  ;  i_leaf < (int)leaves.size()  ;  i_leaf++  ){
			OctreeCell* leaf;   
      leaf = leaves[ i_leaf ];
			double min_coord[ 3 ];
			double max_coord[ 3 ];
			size_t min_key[ 3 ];
			size_t max_key[ 3 ];
   	  leaf->GetKey(  6 , max_key  );
   	  leaf->GetKey(  0 , min_key  );
      for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
        min_coord[ i_dim ] = (double)min_key[ i_dim ] / (double)(1<<ROOT_LEVEL);
        max_coord[ i_dim ] = (double)max_key[ i_dim ] / (double)(1<<ROOT_LEVEL);
      }
			leaf->SetIntersection( Intersects( min_coord , max_coord , center , radius  ) );
		}

    //SENDING LEAVES TO SLAVES
    local_root = leaves[ 0 ];
    for(  size_t i_leaf = 1 ; i_leaf < leaves.size()  ;  i_leaf++  ){
      com->MasterSendCellToSlave( (int)i_leaf , leaves[ i_leaf ] );
    }
  
    //SEND POINTS BELONGING TO EACH CELL
    sort( Points.begin( ), Points.end( ), [ ]( Node& lhs, Node& rhs ){ return lhs._GetRank() < rhs._GetRank(); });
    for(  size_t i_leaf = 1 ; i_leaf < leaves.size()  ;  i_leaf++  ){
      com->MasterSendNodesToSlave( (int)i_leaf , Points );
    }
    
    std::vector<Node> Points_copy;
    Points_copy.push_back( Points[ 0 ] );
    int actual_rank = Points[ 0 ]._GetRank();
    int next_rank;
    int *begin_indexes = new int [com->_GetSize()];
    int flag = 1;
    begin_indexes[ 0 ] = 0;
    for(  size_t i_node = 1  ;  i_node < Points.size()  ;  i_node++  ){
      Points_copy.push_back( Points[ i_node ] );
      next_rank = Points[ i_node ]._GetRank();
      if(  next_rank != actual_rank  ){
        begin_indexes[ flag ] = i_node;
        flag++;
      }
      actual_rank = next_rank;
    }


    //RANK MASTER PROCESS LEAF 0 AND POINTS WITH RANK 0 (ERASING EXTRA INFORMATION)
    int n_nodes = (int)Points.size();
    int i_point = 0;
    while(  i_point < ( n_nodes - ( begin_indexes[ 1 ] - begin_indexes[ 0 ] ) )  ){
      Points.pop_back();
      i_point++;
    }
    leaves.clear();



    //SUBDIVIDING CELLS BASED ON INTERSECTION WITH THE SPHERE
    int i_level = 0;
    while( i_level < level ){
      i_level++;
      //OBTAINING ALL LEAVES
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
      //SUBDIVIDING UNTIL REFINEMENT LEVEL
			OctreeCell* leaf;      
			const int n_leaves = (int)leaves.size();
			for(  int i_leaf = 0  ;  i_leaf < n_leaves  ;  i_leaf++  ){
				leaf = leaves[ i_leaf ];
				if( leaf->GetIntersection() ){
					leaf->SubdivideCell(  );
          //CALCULATING THE CHILDS INTERSECTED BY CIRCLE
          OctreeCell* child;
          for(  size_t i_child = 0  ;  i_child < 8  ;  i_child++  ){
            child = leaf->pGetChild( i_child );
				    size_t min_key[ 3 ];
				    size_t max_key[ 3 ];
				    double max_coord[ 3 ];
				    double min_coord[ 3 ];
        	  child->GetKey(  6 , max_key  );
        	  child->GetKey(  0 , min_key  );
            for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
              max_coord[ i_dim ] = (double)max_key[ i_dim ] / (double)(1<<ROOT_LEVEL);
              min_coord[ i_dim ] = (double)min_key[ i_dim ] / (double)(1<<ROOT_LEVEL);
            }
      			child->SetIntersection( Intersects( min_coord , max_coord , center , radius  ) );
          }
				}
			}
			leaves.clear();
    }

    //OBTAINING ALL LEAVES AND THE NODES THAT INTERSECTS THE SPHERE
    std::vector<int> indexes;
    for(  size_t i_node = 0  ;  i_node < Points.size()  ;  i_node++  ){
      key_type keys[ 3 ];
      for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  )
        keys[ i_dim ] = Points[ i_node ]._GetKey( i_dim );

      OctreeCell* cell = local_root;
      for(  size_t i = 0; i < (ROOT_LEVEL-1) ; i++  ){
        if( cell->IsLeaf() ) {
          break;
        }
        cell = cell->pGetChild( keys[ 0 ] , keys[ 1 ] , keys[ 2 ] );
      }
      if(  cell->GetIntersection()  ){

        double dist = 0.0;
        for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  )
          dist += pow( Points[ i_node ]._GetCoord( i_dim ) - center[ i_dim ] , 2 );  
        
        dist = pow( dist , 0.5 );
        if(  dist <= radius  ){
          indexes.push_back( (int)i_node );
        }
      }
    }
    std::cout<<"Master encontro: "<<indexes.size()<<std::endl;
    delete[] begin_indexes;
  }else{//SLAVES

    std::vector<Node> Points;
    //RECEIVING PROCCESS INFO
    local_root=com->SlaveReceiveCellFromMaster( );
    Points = com->SlaveReceiveNodesFromMaster( );

    //SUBDIVIDING AND CALCULATING INTERSECTIONS
    OctreeCell_vector leaves;
    int i_level = 0;
    while( i_level < level ){
      i_level++;
      //OBTAINING ALL LEAVES
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

      //SUBDIVIDING UNTIL REFINEMENT LEVEL
			OctreeCell* leaf;      
			const int n_leaves = (int)leaves.size();
			for (  int i_leaf = 0  ;  i_leaf < n_leaves  ;  i_leaf++  ){
				leaf = leaves[ i_leaf ];
				if( leaf->GetIntersection() ){
					leaf->SubdivideCell(  );
          //CALCULATING THE CHILDS INTERSECTED BY CIRCLE
          OctreeCell* child;
          for(  size_t i_child = 0  ;  i_child < 8  ;  i_child++  ){
            child = leaf->pGetChild( i_child );
				    size_t min_key[ 3 ];
				    size_t max_key[ 3 ];
				    double min_coord[ 3 ];
				    double max_coord[ 3 ];
        	  child->GetKey(  0 , min_key  );
        	  child->GetKey(  6 , max_key  );
            for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
              max_coord[ i_dim ] = (double)max_key[ i_dim ] / (double)(1<<ROOT_LEVEL);
              min_coord[ i_dim ] = (double)min_key[ i_dim ] / (double)(1<<ROOT_LEVEL);
            }
      			child->SetIntersection( Intersects( min_coord , max_coord , center , radius ) );
          }
				}
			}
			leaves.clear();
    }
    
    //OBTAINING ALL LEAVES AND THE NODES THAT INTERSECTS THE SPHERE    
    std::vector<int> indexes;
    for(  size_t i_node = 0  ;  i_node < Points.size()  ;  i_node++  ){
      key_type keys[ 3 ];
      for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  )
        keys[ i_dim ] = Points[ i_node ]._GetKey( i_dim );

      OctreeCell* cell = local_root;
      for(  size_t i = 0; i < (ROOT_LEVEL-1) ; i++  ){
        if(cell->IsLeaf()) {
          break;
        }
        cell = cell->pGetChild( keys[ 0 ] , keys[ 1 ] , keys[ 2 ] );
      }
      if(  cell->GetIntersection()  ){

        double dist = 0.0;
        for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  )
          dist += pow( Points[ i_node ]._GetCoord( i_dim ) - center[ i_dim ] , 2 );  
        
        dist = pow( dist , 0.5 );
        if(  dist <= radius  ){
          indexes.push_back( (int)i_node );
        }
      }
    }
    std::cout<<com->_GetRank()<<" Esclavo encontro: "<<indexes.size()<<std::endl;
  }
	MPI_Finalize();
}






