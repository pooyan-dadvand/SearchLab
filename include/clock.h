// FemLab -*- C++ -*- Created: galindo Sep-1993, Modified: galindo Sep-1993
#ifndef Clock_h
#define Clock_h

#ifdef _WIN32
#include <time.h>
#include <sys/types.h>
#include <sys/timeb.h>
#else
#include <unistd.h>
#include <sys/times.h>
#include <sys/timeb.h>
#include <time.h> // CLK_TCK
#  ifndef CLK_TCK
#   define CLK_TCK      CLOCKS_PER_SEC
#  endif
#endif

//class ofstream; // defined elsewhere
//class ostream;

#include <stdio.h>

typedef enum  {
  CLOCK_NICE_TEXT_FORMAT = 0,
  CLOCK_TAB_SEPARATED_VALUES = 1
} t_ClockFormatEnum;

// **********************************
class TimeTable // the CPU time table
// **********************************
{
 private:

  friend class Clock;

  static TimeTable  *current; // the current link into the list of time tables
  static int        level;
  //static ofstream*  report;
  static FILE       *report;
  static t_ClockFormatEnum m_format;

  TimeTable   *father, *brother, *firstChild, *lastChild;

  char        *myName;
  clock_t     measuredTime;
  double      childrenTime;
  long long int Uses;

  //RAMSAN adding
  static char VariableName[1024];
  static int VariableValue;

              TimeTable              (const char* name, TimeTable* father);
             ~TimeTable              ();

  void        printName              ( FILE *os) const;

 public:

  static void printCurrentStackNames ( FILE *os) ;  //falta const
  static void EnterVariableNameAndValue(char* name,int value);
};

// *********************************************
class Clock // records a time into the TimeTable
// *********************************************
{
 private:
  TimeTable *tablesFather;
  clock_t   clocksBirthday;

  void openTimeFile   (const char* filename) const;

 public:
              Clock   (const char* name);
  // if there is no file, use filename else use entryname
              Clock   (const char* filename,const char* entryname);
             ~Clock   ();

  void EnterVariableNameAndValue(char* name,int value)
       {TimeTable::EnterVariableNameAndValue(name,value);}
  void SetOutputFormat( t_ClockFormatEnum new_format) {
    TimeTable::m_format = new_format;
  }
};

class ConditionalClock {
public:
  ConditionalClock( bool condition, const char *name ) : _ptr_clk( NULL ) {
    if ( condition )
      _ptr_clk = new Clock( name );
  }
  ConditionalClock( bool condition, const char *filename, const char *name ) : _ptr_clk( NULL ) {
    if ( condition )
      _ptr_clk = new Clock( filename, name );
  }
  ~ConditionalClock() {
    if ( _ptr_clk )
      delete _ptr_clk;
  }

private:
  Clock *_ptr_clk;
};

inline clock_t currentTime() {
// #ifndef _WIN32
//   struct tms time;
//   times(&time);
//   // return time.tms_utime; // CPU time
//   // return time.tms_stime; // sys time
//   return time.tms_utime + time.tms_stime; // CPU + sys time
// #else
   return clock();
// #endif
}

// inline time_t currentWallTime() {
//   return time();
// }
inline int currentWallTime( struct timeb *tp) {
#ifndef _WIN32
  return ftime( tp);
#else
  // MS Wind de los co...... que no respeta los estandares.
  ftime( tp);
  return 0;
#endif
}

typedef enum {
  CRONO_SYS_TIME = 0, CRONO_WALL_TIME = 1
} _t_crono_measure_type;

class Crono {
public:
  void ini( void) {
    if ( _time_type == CRONO_SYS_TIME)
      _clock_inicio = currentTime();
    else
      currentWallTime( &_time_inicio);
  }
  float fin( void) {
    float res = 0.0;
    if ( _time_type == CRONO_SYS_TIME)
      res = ( float)( currentTime() - _clock_inicio) / ( float)CLK_TCK;
    else {
      struct timeb tmp;
      currentWallTime( &tmp);
      time_t dif_sec = tmp.time - _time_inicio.time;
      int dif_msec = tmp.millitm - _time_inicio.millitm;
      res = ( float)( ( double)( dif_sec * 1000 + dif_msec) / 1000.0);
    }
    ini();
    return res;
  }
  Crono( _t_crono_measure_type tt = CRONO_SYS_TIME): _time_type( tt) {
    ini();
  };

private:
  _t_crono_measure_type _time_type;
  clock_t _clock_inicio;
  struct timeb _time_inicio;
};

#endif // #ifndef Clock_h
