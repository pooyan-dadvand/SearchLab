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
#include "geometric_operations.h"
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
  int n_neighbours;
  
public:
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

  double _GetCoord( int i_pos ){
    return coords[ i_pos ];
  }

  int _GetRank(){
    return rank;
  }

  key_type _GetKey( int i_pos ){
    return keys[ i_pos ];
  }

  int _GetNNeighbours(){
    return n_neighbours;
  }
  
  void SumNeighbour(  ){
    n_neighbours++;
  }

  void Intersects( double radius , Node* neighbour_node ){
    double neighbour_center[ 3 ];
    neighbour_center[ 0 ] = neighbour_node->_GetCoord( 0 );
    neighbour_center[ 1 ] = neighbour_node->_GetCoord( 1 );
    neighbour_center[ 2 ] = neighbour_node->_GetCoord( 2 );
    if(  Distance( coords , neighbour_center ) <= radius  ){
      n_neighbours++;
      neighbour_node->SumNeighbour();
    }
  }

  void Intersects( double radius , double* neighbour_coords ){
    if(  Distance( coords , neighbour_coords ) <= radius  ){
      n_neighbours++;
    }
  }

};

void ReadPointsInfo( std::string filename , int refinement_level , std::vector<Node*>& Points ){
  std::ifstream myfile;
  //OPENING POINTS DATA FILE
  myfile.open( filename );
	if(!myfile) {
		std::cout << "******ERROR: Can not open points file******" << std::endl;
		exit( 0 );
	}
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
    Node* aux = new Node( coords[ 0 ] , coords[ 1 ] , coords[ 2 ] , rank , keys[ 0 ]+flag[ 0 ] , keys[ 1 ]+flag[ 1 ] , keys[ 2 ]+flag[ 2 ] );
    Points.push_back( aux );
    
  }  
  //CLOSING POINTS DATA FILE
  myfile.close();
}

struct DataInsideOctreeCell
{
  std::vector<Node*> Points;
  std::vector<int> NL_indexes;
public:
  DataInsideOctreeCell( ){

  }
  ~DataInsideOctreeCell(){

  }
  size_t _GetNNodes(){
    return Points.size();
  }

  size_t _GetNNeighbours(){
    return NL_indexes.size();
  }

  int GetIndex( int i_pos ){
    return NL_indexes[ i_pos ];
  }

  void SetNode( Node* node ){
    Points.push_back( node );
  }

  Node* GetNode( size_t i_node ){
    return Points[ i_node ];
  }

  void SetNeighbour( int rank ){
    bool flag = false;
    for(  size_t i_pos = 0  ;  i_pos < NL_indexes.size()  ;  i_pos++   ){
      if(  NL_indexes[ i_pos ] == rank  ){
        flag = true;
        break;
      }
    }
    if(  !flag  ){
      NL_indexes.push_back( rank );
      sort(NL_indexes.begin(), NL_indexes.end(), [](const int a, const int b) {return a > b; });
    }
  }

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

  void CountNeighbourRankNeighbours( int rank , double radius , std::vector<std::vector<double>>& Recv_coords ){
    size_t n_nodes = Points.size();
    for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
      Node* node = Points[ i_node ];
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


class Comunicator{
	int mpi_rank;
	int mpi_size;
  int g_level;
  int **connectivities;


	public:
		Comunicator(int argc,char **argv){
			MPI_Init( &argc , &argv );
			MPI_Comm_rank( MPI_COMM_WORLD , &mpi_rank );
			MPI_Comm_size( MPI_COMM_WORLD , &mpi_size );
      g_level = (int)(log( mpi_size )/log( 8 ));
		}

		~Comunicator(){
      MPI_Finalize();
      for(  int i_ren = 0  ;  i_ren < mpi_size  ;  i_ren++  ){
        delete[] connectivities[ i_ren ];
      }
      delete[] connectivities;
    } 

		int  _GetRank(){
      return mpi_rank;
    }

		int  _GetSize(){
      return mpi_size;
    }

		int  _GetLevel(){
      return g_level;
    }

		void _Synchronize(){
      MPI_Barrier( MPI_COMM_WORLD );
    }

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

    void CreateConnectivitiesMatrix( OctreeCell_vector& leaves ){
      connectivities = new int*[ mpi_size ];
      for(  int i_ren = 0  ;  i_ren < mpi_size  ;  i_ren++  ){
        connectivities[ i_ren ] = new int[ mpi_size ];
      }
      for(  int i_ren = 0  ;  i_ren < mpi_size  ;  i_ren++  ){
          for(  int i_col = 0  ;  i_col < mpi_size  ;  i_col++){
              connectivities[ i_ren ][ i_col ] = 0;
          }
      }
      for(  size_t i_leaf = 0  ;  i_leaf < leaves.size()  ;  i_leaf++  ){
        OctreeCell* leaf = leaves[ i_leaf ];
        for(  size_t i_dir = 0  ;  i_dir < 18  ;  i_dir++  ){
          key_type neighbour_key[ 3 ];
          if(  leaf->GetNeighbourKey( i_dir , neighbour_key )  ){
            int rank = GetRankFromKey( neighbour_key );
            if(  rank != (int)i_leaf  ){
              connectivities[ i_leaf ][ rank ] = 1;
            }
          }
        }
      }
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

    void MasterSendNodesToSlave( int rank , std::vector<Node*>& Points ){
      int n_nodes = 0;
      int begin_position = -1;
      int flag = 0;
      for(  size_t i_node = 0  ;  i_node < Points.size()  ;  i_node++  ){
        if(  Points[ i_node ]->_GetRank() == rank  ){
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
      if(  n_nodes > 0  ){
        void *buff = malloc(52*n_nodes);
        int position = 0;
        for(  int i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
          double coord[ 3 ];
          key_type keys[ 3 ];
          for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
            coord[ i_dim ] = Points[ begin_position + i_node ]->_GetCoord( i_dim ); 
            keys[ i_dim ] = Points[ begin_position + i_node ]->_GetKey( i_dim ); 
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
    }

    void SlaveReceiveNodesFromMaster( std::vector<Node*>& Points ){
      int n_nodes;
      MPI_Recv( &n_nodes , 4 , MPI_BYTE , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
      if(  n_nodes > 0  ){
        void* buff = malloc( 52 * n_nodes );
        MPI_Recv( buff , 52 * n_nodes , MPI_BYTE , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
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
          Node* aux = new Node( coord[ 0 ] , coord[ 1 ] , coord[ 2 ] , rank , keys[ 0 ] , keys[ 1 ] , keys[ 2 ] );
          Points.push_back( aux ); 
        }
        free( buff );
      }
    }

    void AssignBoundaryInformation( OctreeCell* local_root ){
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
          if(  cell->pGetData()  ){
            cell->DeleteData();
          }
        }
      }
    }

    void SetLocalNeighbours( OctreeCell* local_root , double radius ){
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
          data->CountLocalNeighbours( radius );
          for(  size_t i_dir = 0  ;  i_dir < 18  ;  i_dir++  ){
            key_type neighbour_key[ 3 ];
            if(  cell->GetNeighbourKey( i_dir , neighbour_key )  ){
              int rank = GetRankFromKey( neighbour_key );
              if(  rank == mpi_rank   ){
                OctreeCell* neighbour_cell = local_root;
                for(  size_t i_level = 0  ;  i_level < ROOT_LEVEL  ;  i_level++  ){
                    if(  neighbour_cell->IsLeaf()  ) {
                        break;
                    }
                    neighbour_cell = neighbour_cell->pGetChild( neighbour_key[ 0 ] , neighbour_key[ 1 ] , neighbour_key[ 2 ] );
                }
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

    void GetCoordsToSend( OctreeCell_vector& leaves , std::vector<std::vector<double>>& coords ){
      size_t n_leaves = leaves.size();
      for(  size_t i_leaf = 0  ;  i_leaf < n_leaves  ;  i_leaf++  ){
        OctreeCell* cell = leaves[ i_leaf ];
        if(  cell->GetIntersection()  ){
          DataInsideOctreeCell* data = cell->pGetData();
          size_t n_neighbours = data->_GetNNeighbours();
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

void AsignIntersectionOnCells( OctreeCell_vector& leaves , std::vector<Node*>& Points ){
  size_t n_points = Points.size();
  int* flag = new int [ leaves.size() ];
  size_t cont = 0;
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
		  OctreeCell* leaf;
		  size_t n_leaves = leaves.size();
		  for(  size_t i_leaf = 0  ;  i_leaf < n_leaves  ;  i_leaf++  ){
			  leaf = leaves[ i_leaf ];
			  if( leaf->GetIntersection() ){
				  leaf->SubdivideCell(  );
          DataInsideOctreeCell** data_children = new DataInsideOctreeCell*[ 8 ];
          for(  size_t i_child = 0  ;  i_child < 8  ;  i_child++  ){
            data_children[ i_child ] = new DataInsideOctreeCell;
            OctreeCell* child = leaf->pGetChild( i_child );
            DataInsideOctreeCell** data_child = child->pGetDataPointer();
            (*data_child) = data_children[ i_child ]; 
          }

          //ASSIGNING DATA TO CHILDS
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

void RunMultiplePointMPISearchOctree( std::string filename , double radius , int arg , char* argv[] ){
  
  //BEGINING MPI PARALLEL REGION
	Comunicator *com=new Comunicator( arg , argv );
  //THIS IS THE CELL POINTER THAT WILL BE PROCESSED ON EACH MPI PARTITION
  OctreeCell* local_root = NULL;

  //CALCULATING THE REFINEMENT LEVEL IN THE LOCAL OCTREE, THE CELL SIZE HAVE TO BE THE 
  //SMALLEST POSSIBLE BUT IT HAVE NOT TO BE LESS THAN THE RADIUS
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
  //IN THIS CASE IS SUBSTRACTED THE REFINEMENT THAT WILL BE PERFORMED BY THE MASTER
  level = level - com->_GetLevel();


  if( com->_GetRank() == 0 ){

    //MASTER PROCESS
    double t0, t1 , t2; 
    t0 = MPI_Wtime();
    t1 = MPI_Wtime();
    t2 = MPI_Wtime() - t0;
    fprintf(stdout, "[%16lf]-[%16lf]   **MASTER-BEGINNING MULTIPLE POINTS SEARCH\n\n\n", t2 , MPI_Wtime() - t1 );
    
    //CREATING OCTREE AND REFINING UNTIL HAVE A CELL FOR EACH MPI PROCESS
    t1 = MPI_Wtime();
    double cell_size = 1.0 / pow( 2 , com->_GetLevel() );
	  OctreeDriver octree_driver;
	  octree_driver.RefineWithUniformSizeNormalized( cell_size );
    t2 = MPI_Wtime() - t0;
		fprintf(stdout, "[%16lf]-[%16lf]   **Refinement             : \n", t2 , MPI_Wtime() - t1 );

    //OBTAINING ALL LEAVES
    t1 = MPI_Wtime();
	  OctreeCell_vector leaves;
	  leaves = octree_driver.CalcAllLeavesVector( );
    com->CreateConnectivitiesMatrix( leaves );
    t2 = MPI_Wtime() - t0;
		fprintf(stdout, "[%16lf]-[%16lf]   **Calculating all leaves : \n", t2 , MPI_Wtime() - t1 );

    //READING POINTS INFO
    t1 = MPI_Wtime();
    std::vector<Node*> Points;
    ReadPointsInfo( filename , com->_GetLevel() , Points );
    t2 = MPI_Wtime() - t0;
		fprintf(stdout, "[%16lf]-[%16lf]   **Read points information: \n", t2 , MPI_Wtime() - t1 );

    //DETERMINING LEAVES INTERSECTED BY POINTS
    t1 = MPI_Wtime();
    AsignIntersectionOnCells( leaves , Points );
    t2 = MPI_Wtime() - t0;
		fprintf(stdout, "[%16lf]-[%16lf]   **Setting intersections  : \n", t2 , MPI_Wtime() - t1 );

    //SENDING LEAVES TO SLAVES
    t1 = MPI_Wtime();
    local_root = leaves[ 0 ];
    for(  size_t i_leaf = 1 ; i_leaf < leaves.size()  ;  i_leaf++  ){
      com->MasterSendCellToSlave( (int)i_leaf , leaves[ i_leaf ] );
    }
    t2 = MPI_Wtime() - t0;
		fprintf(stdout, "[%16lf]-[%16lf]   **Sending leaves         : \n", t2 , MPI_Wtime() - t1 );

    //SEND POINTS BELONGING TO EACH CELL
    t1 = MPI_Wtime();
    sort( Points.begin( ), Points.end( ), [ ]( Node* lhs, Node* rhs ){ return lhs->_GetRank() < rhs->_GetRank(); });
    for(  size_t i_leaf = 1 ; i_leaf < leaves.size()  ;  i_leaf++  ){
      com->MasterSendNodesToSlave( (int)i_leaf , Points );
    }
    CleanLeavesAndPoints( leaves , Points );
    t2 = MPI_Wtime() - t0;
		fprintf(stdout, "[%16lf]-[%16lf]   **Sending points         : \n", t2 , MPI_Wtime() - t1 );


    //SUBDIVIDING LOCAL ROOT 
    t1 = MPI_Wtime(); 
    RefineCellsIntersected( local_root , Points , level );
    t2 = MPI_Wtime() - t0;
		fprintf(stdout, "[%16lf]-[%16lf]   **Refining octree        : \n", t2 , MPI_Wtime() - t1 );

    //OBTAINING BOUNDARY INFORMATION  
    t1 = MPI_Wtime();
    com->AssignBoundaryInformation( local_root );
    t2 = MPI_Wtime() - t0;
		fprintf(stdout, "[%16lf]-[%16lf]   **Setting boundary info  : \n", t2 , MPI_Wtime() - t1 );

    //NEIGHBOURS SEARCH
    t1 = MPI_Wtime();
    com->NeighboursSearch( local_root , radius );
    t2 = MPI_Wtime() - t0;
		fprintf(stdout, "[%16lf]-[%16lf]   **Neighbours search      : \n", t2 , MPI_Wtime() - t1 );  

    //PRINTING TOTAL TIME
    fprintf(stdout, "                                        **_______________________:\n");
    t2 = MPI_Wtime() - t0;
		fprintf(stdout, "[%16lf]                      **TOTAL TIME             : \n", t2 );
  }else{
    //SLAVES
    std::vector<Node*> Points;

    //RECEIVING INFORMATION FROM MASTER
    local_root=com->SlaveReceiveCellFromMaster( );
    com->SlaveReceiveNodesFromMaster( Points );

    //SUBDIVIDING LOCAL ROOT
    RefineCellsIntersected( local_root , Points , level );

    //SETTING BOUNDARY INFORMATION
    com->AssignBoundaryInformation( local_root );

    //NEIGHBOURS SEARCH
    com->NeighboursSearch( local_root , radius );

  }
	MPI_Finalize();
}






