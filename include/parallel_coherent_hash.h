/* -*- c++ -*- */
#pragma once

#include <array>
#include <random>
#include <vector>

#include <iostream> // for PrintStatistics()
#include <locale> // for thousand separators

#include "primes.h"
#include "sl_exception.h"

template < typename TDataType, typename t_KeyType = std::uint64_t >
class ParallelCoherentHash {
private:
  static constexpr int MAXIMUM_AGE = 63;
  static constexpr double LOAD_FACTOR = 0.9;

public:
  // using t_KeyType = std::uint64_t;

private:
  // internal class definitions
  t_KeyType getRandomOffset() {
    std::random_device rd; // rd() returns a random 'unsigned int'
    constexpr size_t size_random_value = sizeof( unsigned int);
    size_t num_random = ( sizeof( t_KeyType) + sizeof( t_KeyType) - 1) / size_random_value;
    t_KeyType value = 0;
    for ( size_t i = 0; i < num_random; i++) {
      value = ( value << ( i * size_random_value)) ^ ( t_KeyType)rd();
    }
    if ( value < 0) // may be t_KeyType is signed, and we want positive offsets
      value = -value;
    return value;
  } 
  
  class KeyAgeTuple {
  public:
    KeyAgeTuple() {}
    KeyAgeTuple( const KeyAgeTuple &kat): m_Key( kat.m_Key), m_Age( kat.m_Age) {}
    KeyAgeTuple( const t_KeyType &key_in, const char age_in = 0): m_Key( key_in), m_Age( age_in) {}
    KeyAgeTuple &operator=( const KeyAgeTuple &kat) { m_Key = kat.m_Key; m_Age = kat.m_Age; return *this;}
      t_KeyType getKey() const { return m_Key;}
    char getAge() const { return m_Age;}
    void setKey( const t_KeyType &key_in) { m_Key = key_in;}
    void setAge( const char &age_in) { m_Age = age_in;}
  private:
    t_KeyType m_Key;
    char m_Age;
  };
  
  class HashEntry {
    char m_MaximumAge; // if m_MaximumAge == 0, then it's empty
    KeyAgeTuple m_KeyAge;
    TDataType m_Data;
  public:
    HashEntry() : m_MaximumAge( 0 ), m_KeyAge( 0 ), m_Data() {}
    HashEntry( const t_KeyType &key_in, const TDataType &data_in )
      : m_MaximumAge( 0 ), m_KeyAge( key_in ), m_Data( data_in ) {}
    HashEntry( HashEntry const &he )
      : m_MaximumAge( he.m_MaximumAge), m_KeyAge( he.m_KeyAge), m_Data( he.m_Data) {}
    HashEntry &operator=( const HashEntry &he) {
      m_MaximumAge = he.m_MaximumAge;
      m_KeyAge = he.m_KeyAge;
      m_Data = he.m_Data;
      return *this;
    }
    void copyKeyAgeData(  const HashEntry &he) {
      // m_MaximumAge = he.m_MaximumAge; // we want to maintain the maximum age!
      m_KeyAge = he.m_KeyAge;
      m_Data = he.m_Data;
    }
    void setMaximumAge( const char max_age_in ) { m_MaximumAge = max_age_in; }
    char getMaximumAge() const { return m_MaximumAge; }

    t_KeyType getKey() const { return m_KeyAge.getKey(); }
    void setKey( const t_KeyType &key_in ) { m_KeyAge.setKey( key_in); }
    char getAge() const { return m_KeyAge.getAge(); }
    void setAge( const char age_in ) { m_KeyAge.setAge( age_in); }
    const TDataType &getData() const { return m_Data; }
    TDataType &getDataRef() { return m_Data; }
    void setData( const TDataType &data_in ) { m_Data = data_in; }
  };

  // ParallelCoherentHash methods and members
public:
  ParallelCoherentHash(): m_Size( 0), m_TableMaximumAge( 0) { FillRandomOffsets(); }
  void resize( std::size_t num_entries, const TDataType &initial_value);
  void Insert( const t_KeyType &key_in, const TDataType &data_in, bool &found_or_inserted_out );
  
  // const TDataType &operator[]( const t_KeyType &key_in) const {
  const TDataType &getData( const t_KeyType &key_in) const {
    bool found = false;
    const TDataType &ret = getData( key_in, found);
    if ( !found) {
      // throw SLexception( "ParallelCoherentHash::operator[] : couldn't find key = " + std::to_string( key_in));
      return m_InitialValue;
    }
    return ret;
  }
  // TDataType &operator[]( const t_KeyType &key_in) {
  TDataType &getDataRef( const t_KeyType &key_in) {
    bool found_or_inserted = false;
    TDataType &ret = getDataRef( key_in, found_or_inserted);
    if ( !found_or_inserted) {
      throw SLexception( "ParallelCoherentHash::operator[] : couldn't find or insert key = " + std::to_string( key_in));
    }
    return ret;
  }

  void PrintStatistics() const;

  // for statistic purposes:
  std::size_t getRawHashTableSize() const { return m_HashTable.size();}
  TDataType getRawEntryData( std::size_t idx) const { return m_HashTable[ idx].getData();}
  TDataType &getRawEntryDataRef( std::size_t idx) { return m_HashTable[ idx].getDataRef();}
  bool getRawEntryUsed( std::size_t idx) const { return m_HashTable[ idx].getAge() != 0;}
  
private:
  const TDataType &getData( const t_KeyType &key_in, bool &found_out) const;
  TDataType &getDataRef( const t_KeyType &key_in, bool &found_or_inserted_out);
  TDataType &findDataRef( const t_KeyType &key_in, bool &found_out);
  void FillRandomOffsets() {
    for ( auto &offset : m_HashOffsets ) {
      offset = this->getRandomOffset();
    }
    m_HashOffsets[ 0] = 0; // First offset = 0
    // for ( size_t i = 0; i < m_HashOffsets.size(); i++) {
    //   std::cout << "Random offset " << i << " = " << m_HashOffsets[ i] << std::endl;
    // }
  }
  std::size_t getProperSize( std::size_t num_entries); // returns next "prime" above num_entries/load_factor
  std::size_t getHash( const t_KeyType &key_in, const char age_in) const {
    // doing modulus with std::size_t we don't care if t_KeyType is signed and key_in may be < 0
    return ( std::size_t)( key_in + m_HashOffsets[ age_in]) % m_Size;
  }
  std::size_t countUsedEntries() const {
    std::size_t count = 0;
    for ( const auto &entry : m_HashTable ) {
      count += ( entry.getAge() != 0);
    }
    return count;
  }
  
private:
  std::size_t m_Size;
  std::vector< HashEntry > m_HashTable;
  std::array< t_KeyType, MAXIMUM_AGE > m_HashOffsets;
  TDataType m_InitialValue;
  char m_TableMaximumAge;
};

// ParallelCoherentHash member methods:

template < typename TDataType, typename t_KeyType>
inline std::size_t ParallelCoherentHash< TDataType, t_KeyType>::getProperSize( std::size_t num_entries) {
  std::size_t desired_size = ( std::size_t)( 0.5 + ( double)num_entries / LOAD_FACTOR);
  std::size_t ret = 0;
  if ( ( std::size_t)get_last_prime() > desired_size) {
    ret = ( std::size_t)get_ceil_prime( ( unsigned int)desired_size); // get_last_prime is an unsigned int
  } else {
    std::size_t prim = ( std::size_t)next_highest_square_prime( ( unsigned long long int)desired_size); // we only know primes up to 3M - 1,
    ret = prim * prim;
  }
  return ret; // there is no minimum size ...
}

template < typename TDataType, typename t_KeyType>
inline void ParallelCoherentHash< TDataType, t_KeyType>::resize( std::size_t num_entries, const TDataType &initial_value) {
  m_InitialValue = initial_value;
  if ( m_Size == 0) {
    m_Size = this->getProperSize( num_entries);
    HashEntry val( ( t_KeyType)0, initial_value);
    m_HashTable.resize( m_Size, val);
    m_TableMaximumAge = 0;
  } else {
    throw SLexception( "ParallelCoherentHash::resize : Resizing a non empty hash table is not supported yet\n" );
  }
}

template < typename TDataType, typename t_KeyType>
inline void ParallelCoherentHash< TDataType, t_KeyType>::Insert( const t_KeyType &key_in, const TDataType &data_in,
								 bool &found_or_inserted_out) {
  found_or_inserted_out = false;
  char age = 1;
  HashEntry entry_to_insert( key_in, data_in);
  // if ( countUsedEntries() > m_Size - 1) {
  //   std::cout << "ParallelCoherentHash::Insert TABLE FULL inserting key = " << key_in << std::endl;
  //   std::cout << "ParallelCoherentHash::Insert TABLE of size " << m_Size << " has " << countUsedEntries() << " entries " << std::endl;
  // }
  while ( age < MAXIMUM_AGE) {
    entry_to_insert.setAge( age);
    std::size_t idx = getHash( entry_to_insert.getKey(), ( char)( age - 1));
    HashEntry old_entry( m_HashTable[ idx]);
    // if entry_to_insert is older than stored one, evict it!
    if ( entry_to_insert.getAge() > old_entry.getAge()) {
      // we want to maintain the maximum age of the entry
      m_HashTable[ idx].copyKeyAgeData( entry_to_insert);
    }
    if ( entry_to_insert.getAge() > old_entry.getAge()) { // eviction occured
      // update MaxAge with current age if greated than stored on first insert
      std::size_t age0entryIdx = getHash( entry_to_insert.getKey(), 0);
      // atomicMAX
      char old_age0entry = m_HashTable[ age0entryIdx].getMaximumAge();
      if ( entry_to_insert.getAge() > old_age0entry) {
	m_HashTable[ age0entryIdx].setMaximumAge( entry_to_insert.getAge());
	// global statistics:
	if ( entry_to_insert.getAge() > m_TableMaximumAge)
	  m_TableMaximumAge = entry_to_insert.getAge();
      }
      if ( old_entry.getAge() > 0) {
	// evicted cell was not empty, so process evicted key
	entry_to_insert.copyKeyAgeData( old_entry);
	age = old_entry.getAge();
	age++; // next try for evicted element
      } else {
	// the evicted cell was empty, so quit
	break;
      }
    } else {
      // no eviction occurred =
      // entry_to_insert still younger than key in table, so increase age and try again
      age++;
    }
  } // while
  if ( age >= MAXIMUM_AGE) {
    throw SLexception( "ParallelCoherentHash::Insert couldn't insert key = " + std::to_string( key_in));
  }
  found_or_inserted_out = ( age != MAXIMUM_AGE);
} // ParallelCoherentHash::Insert

template < typename TDataType, typename t_KeyType>
inline TDataType &ParallelCoherentHash< TDataType, t_KeyType>::findDataRef( const t_KeyType &key_in, bool &found_out) {
  // TDataType &ret = m_InitialValue;
  char age = 0;
  std::size_t idx= this->getHash( key_in, age);
  char max_age = m_HashTable[ idx].getMaximumAge(); // maximum number of tries
  t_KeyType key_found = m_HashTable[ idx].getKey();
  if ( m_HashTable[ idx].getAge()) { // cell is not empty
    if ( key_found == key_in) {
      // key found in the first attempt
      // ret = m_HashTable[ idx].getDataRef();
      found_out = true;
      return m_HashTable[ idx].getDataRef();
    } else {
      // cell is not empty and key_in is not there
      char iage = 1;
      for ( iage = 1; iage < max_age; iage++) {
	idx = this->getHash( key_in, iage);
	key_found = m_HashTable[ idx].getKey();
	if ( key_found == key_in) {
	  // ret = m_HashTable[ idx].getDataRef();
	  found_out = true;
	  return m_HashTable[ idx].getDataRef();
	  // break;
	}
      }
      if ( iage == max_age) { // not found
	found_out = false;
      }
    } // if key_in == key_found & else
  } else {
    // cell is empty, so key not found
    found_out = false;
  }
  return m_InitialValue; // somewhat dangerous !!!
}

template < typename TDataType, typename t_KeyType>
inline const TDataType &ParallelCoherentHash< TDataType, t_KeyType>::getData( const t_KeyType &key_in, bool &found_out) const {
  return ( const TDataType &)( ( ParallelCoherentHash< TDataType, t_KeyType> *)this)->findDataRef( key_in, found_out);
}

template < typename TDataType, typename t_KeyType>
inline TDataType &ParallelCoherentHash< TDataType, t_KeyType>::getDataRef( const t_KeyType &key_in, bool &found_or_inserted_out) {
  found_or_inserted_out = false;
  {
    TDataType &ret = findDataRef( key_in, found_or_inserted_out);
    if ( found_or_inserted_out) return ret;
  }
  Insert( key_in, m_InitialValue, found_or_inserted_out);
  if ( found_or_inserted_out) {
    return findDataRef( key_in, found_or_inserted_out);
  }
  found_or_inserted_out = false;
  return m_InitialValue; // somewhat dangerous !!!
}

template < typename TDataType, typename t_KeyType>
inline void ParallelCoherentHash< TDataType, t_KeyType>::PrintStatistics() const {
  std::cout << "=== ParallelCoherentHash statistics === \n";
  std::locale prev_loc = std::cout.getloc();
  std::cout.imbue( std::locale( "")); // for thousand separators ...
  std::cout << "sizeof t_KeyType = " << sizeof( t_KeyType) << std::endl;
  std::cout << "sizeof TDataType = " << sizeof( TDataType) << std::endl;
  std::cout << "sizeof HashEntry = " << sizeof( HashEntry) << std::endl;
  std::size_t sizeHashTable = m_Size * sizeof( HashEntry);
  std::size_t sizeOffsetTable = m_HashOffsets.size() * sizeof( t_KeyType);
  std::cout << "Total number of elements in hash table = " << m_Size << " = " << sizeHashTable << " bytes" << std::endl;
  std::size_t count = 0;
  std::cout << "Hash entry XXXX = (  key, data, age)" << std::endl;
  std::cout.imbue( prev_loc); // restore previous locale, i.e. without thousand separators
  for ( const auto &entry : m_HashTable ) {
    count += ( entry.getAge() != 0);
    // if ( entry.getAge()) {
    //   std::cout << "Hash entry " << count << " = ( " << entry.getKey() << ", " << entry.getData() << ", " << ( int)( entry.getAge()) << ")" << std::endl;
    // }
  }
  std::cout.imbue( std::locale( "")); // for thousand separators ...
  std::cout << " Used number of elements in hash table = " << count << std::endl;
  std::cout << "Configured load factor = " << LOAD_FACTOR << " vs real load factor = " << ( double)count / ( double)m_Size << std::endl;
  std::cout << "Maximum age = " << ( int)m_TableMaximumAge << std::endl;
  std::cout << "Total size in bytes = " << sizeof( *this) + sizeHashTable + sizeOffsetTable << std::endl;
  std::cout.imbue( prev_loc); // restore previous locale, i.e. without thousand separators
  std::cout << "=== End of statistics === \n";
}
