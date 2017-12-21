
// FemLab -*- C++ -*- Created: galindo Sep-1993, Modified: galindo Sep-1993
#include "clock.h"

#include <stdio.h>
#include <stdlib.h>

#if ( defined(__GNUC__) && (__GNUC__>2)) || ( defined(_MSC_VER) && ( _MSC_VER >= 1310))
#include <fstream>
#include <iomanip>
using std::ofstream;
#else
#include <fstream.h>   // report
#include <iomanip.h>   // setw()
#endif
#include <string.h>    // strlen
#ifndef _WIN32
#include <unistd.h>
#include <sys/param.h> // CLK_TCK
#endif
#include <time.h> // CLK_TCK

#define FPRINTF_UTF     fprintf
#define PRINTF_UTF      printf
#define Error           printf
#define FOPEN_UTF       fopen
#define _( a)           ( a)
#define GIDstrcmp       strcmp

#define NAME_LENGTH 44

TimeTable  *TimeTable::current = NULL;
int        TimeTable::level    = 1;
FILE       *TimeTable::report  = NULL;
t_ClockFormatEnum TimeTable::m_format = CLOCK_NICE_TEXT_FORMAT;

//RAMSAN adding
int TimeTable::VariableValue = 0;
char TimeTable::VariableName[ 1024];

// static clock_t currentTime()
// {
// #ifndef _WIN32
//   struct tms time;
//   times(&time);
//   return time.tms_utime; // CPU time
// #else
//   return clock();
// #endif
// }

static void Print( FILE *ofs, int l, const char* n,
                  double t, double p, long long int u,int VariableValue,
                   t_ClockFormatEnum format) // = CLOCK_NICE_TEXT_FORMAT)
{
  if ( format == CLOCK_NICE_TEXT_FORMAT) {
    
    FPRINTF_UTF( ofs, "   ");
    for(int i=2; i<l; i++) FPRINTF_UTF( ofs, " | ");
    if ( l > 1)
      FPRINTF_UTF( ofs, " + ");
    FPRINTF_UTF( ofs, "%s ", n);
    for ( int f = 0; f < (int)(NAME_LENGTH-strlen(n)); f++) FPRINTF_UTF( ofs, ".");
    FPRINTF_UTF( ofs, " %.3g s (%.3g %%)", t, p);
    if(u!=0) FPRINTF_UTF( ofs, " %lld", u);
    else FPRINTF_UTF( ofs, "  ");
    
    //RAMSAN adding
    if(VariableValue && t>0.0){
      //ofs.precision(0);
      FPRINTF_UTF( ofs, " -- %.3g", ( double)VariableValue/t*60.0);
    }
    
    FPRINTF_UTF( ofs, "\n");
    
  } else if ( format == CLOCK_TAB_SEPARATED_VALUES) {
    
    // 1st field: name
    char *name_to_print = new char[ l * 5 + strlen( n)];
    name_to_print[ 0] = '\0';
    for(int i=2; i<l; i++) strcat( name_to_print, " | ");
    if ( l > 1)
      strcat( name_to_print, " + ");
    strcat( name_to_print, n);
    FPRINTF_UTF( ofs, "%s\t", name_to_print);

    // 2nd and 3rd: time(s), time(%)
    FPRINTF_UTF( ofs, "%.3g\t%.3g %%\t", t, p);
    
    // 4th field: # called times
    if(u!=0) FPRINTF_UTF( ofs, "%lld", u);
    else FPRINTF_UTF( ofs, "  ");

    // 5th field: optional variable name
    //RAMSAN adding
    if(VariableValue && t>0.0){
      //ofs.precision(0);
      FPRINTF_UTF( ofs,"\t %.3g", ( double)VariableValue/t*60.0);
    }

    // eol
    FPRINTF_UTF( ofs, "\n");
  }
  
}

// *******************************************
//class TimeTable - the private CPU time table
// *******************************************

TimeTable::TimeTable(const char* n, TimeTable* f) : father(f)
{
  myName = strdup( n);
  brother = firstChild = lastChild = 0;
  measuredTime = 0 ;
  childrenTime = 0.0;
  Uses = 0;
}

TimeTable::~TimeTable()
{
  double myTime = double(measuredTime)/CLK_TCK;
  double percent = 100.0;
  if(father)
    {
      father->childrenTime += myTime;
      double myFathersTime = double(father->measuredTime)/CLK_TCK;
      if(myFathersTime!=0) percent = myTime*100/myFathersTime;
    }

  if ( report) {
    Print(report, level, myName, myTime, percent, Uses,VariableValue, m_format);
  }

  level++;
  if ( firstChild) delete firstChild;
  firstChild = lastChild = 0;

  if(myTime>1.0 && childrenTime!=0.0)
    {
      double otherOperations = myTime-childrenTime;
      percent = otherOperations*100/myTime;
      if(otherOperations>0.1 && percent>0.1) {
	if ( report) {
	  Print(report, level, _("Other operations"), otherOperations, percent, 0,
		VariableValue, m_format);
	}
      }
    }

  level--;
  if ( brother) delete brother;
  brother = 0;
  if ( myName) free( myName);
  myName = NULL;
}

void TimeTable::printName( FILE *os) const
{
  FPRINTF_UTF( os, "\t%s\n", myName);
  if(father) father->printName(os);
}

void TimeTable::printCurrentStackNames( FILE *os) //const
{
  if(TimeTable::current)
    {
      FPRINTF_UTF( os, " Current time table stack:\n\n");
      TimeTable::current->printName(os);
      FPRINTF_UTF( os,"\n");
      if ( report) {
	FPRINTF_UTF( report, "\n\n *** Time Table not produced due to error ***\n\n");
	fclose( report);
      }
      //delete report;
      report = NULL;
    }
}

//RAMSAN adding
void TimeTable::EnterVariableNameAndValue(char* name,int value)
{
  strncpy(VariableName,name,1024);
  VariableName[1023]='\0';
  VariableValue=value;
}

// **********************************************
//class Clock - records a time into the TimeTable
// **********************************************

void Clock::openTimeFile(const char* name) const
{
  int len = ( int)strlen(name);
  char* fname;
  fname = new char[ len + 100];
  strcpy(fname, name);
  strcat( &fname[ len], "_time.txt");

  TimeTable::report = FOPEN_UTF( fname, "w");
  //PRINTF_UTF( "new report %x.\n", TimeTable::report);
#ifdef CLOCK_VERBOSE
  PRINTF_UTF( "Clock file '%s' created", fname);
#endif
  //if( TimeTable::report->fail()) {
  if ( !TimeTable::report || ferror( TimeTable::report)) {
    Error(_("Cannot open output file '%s'"),fname);
    TimeTable::report = NULL;
  }
  delete [] fname;
}

Clock::Clock(const char *name_in) {
  TimeTable *table = NULL;
  if(!TimeTable::current) {
    openTimeFile( name_in);
  }

#if 0
  char name[ 1024];

  // si hay path, le quito el path
  int n = strlen( name_in);
  const char *s = name_in;
  // n - 2 para que al menos haya un caracter
  for ( int i = n - 2; i >= 0; i--) {
    s = &name_in[ i];
    if ( ( *s == '\\') || ( *s == '/')) {
      s++;
      break;
    }
  }
  strcpy( name, s);
#else
  const char *name = &name_in[ strlen( name_in)];
  while ( name > name_in) {
    if ( ( *name == '\\') || ( *name == '/')) {
      name++;
      break;
    }
    name--;
  }
#endif

  tablesFather = TimeTable::current;

  if(tablesFather) table = tablesFather->firstChild;
  while(table && ( GIDstrcmp(table->myName,name))!=0) table = table->brother;
  if(!table)
    {
      table = new TimeTable(name, tablesFather);
      if(tablesFather)
        {
          if(tablesFather->lastChild) tablesFather->lastChild->brother = table;
          else                        tablesFather->firstChild = table;
          tablesFather->lastChild = table;
        }
    }
  table->Uses++;

  TimeTable::current = table;

  clocksBirthday = currentTime();
}

// if there is no file, use filename else use entryname
Clock::Clock(const char* filename,const char* entryname) {
  TimeTable *table = NULL;
  const char* name_in;

  if(!TimeTable::current)
    {
      name_in=filename;
      openTimeFile(name_in);
    }
  else name_in=entryname;

  // as in Clock::Clock( const char *name_in)
  // use only the filename, not the full path
  const char *name = &name_in[ strlen( name_in)];
  while ( name > name_in) {
    if ( ( *name == '\\') || ( *name == '/')) {
      name++;
      break;
    }
    name--;
  }
  
  tablesFather = TimeTable::current;

  if(tablesFather) table = tablesFather->firstChild;
  while(table && strcmp(table->myName,name)!=0) table = table->brother;
  if(!table)
    {
      table = new TimeTable(name, tablesFather);
      if(tablesFather)
        {
          if(tablesFather->lastChild) tablesFather->lastChild->brother = table;
          else                        tablesFather->firstChild = table;
          tablesFather->lastChild = table;
        }
    }
  table->Uses++;

  TimeTable::current = table;

  clocksBirthday = currentTime();
}


Clock::~Clock() {
  clock_t clocksDieday = currentTime();
  
  clock_t clocksLife = clocksDieday - clocksBirthday;

  TimeTable::current->measuredTime += clocksLife;
  
  if(tablesFather) {
    TimeTable::current = tablesFather;
  } else {
    if ( TimeTable::report) {
      if(!TimeTable::VariableValue){
        if ( TimeTable::m_format == CLOCK_NICE_TEXT_FORMAT) {
          fprintf( TimeTable::report, "\n   Time Table: ................................   cputime[s] (percent[%%]) Uses\n\n");
        } else if ( TimeTable::m_format == CLOCK_TAB_SEPARATED_VALUES) {
          // 1st field: name
          // 2nd and 3rd: time(s), time(%)
          // 4th field: # called times
          // 5th field: optional variable name
          fprintf( TimeTable::report, "Time Table\tcputime[s]\tpercent[%%]\tUses\n");
        }
      }
      else{  //RAMSAN adding
        if ( TimeTable::m_format == CLOCK_NICE_TEXT_FORMAT) {
          fprintf( TimeTable::report, "\n Number of %s=%d\n\n", TimeTable::VariableName, TimeTable::VariableValue);
          
          fprintf( TimeTable::report, "\n   Time Table: ................................   "
                                            "cputime[s] (percent[%%]) Uses -- %s/min\n\n", TimeTable::VariableName);
        } else if ( TimeTable::m_format == CLOCK_TAB_SEPARATED_VALUES) {
          // 1st field: name
          // 2nd and 3rd: time(s), time(%)
          // 4th field: # called times
          // 5th field: optional variable name
          fprintf( TimeTable::report, "Number of %s=%d\n", TimeTable::VariableName, TimeTable::VariableValue);
          fprintf( TimeTable::report, "Time Table\tcputime[s]\tpercent[%%]\tUses\t%s/min\n", TimeTable::VariableName);
        }
      }
    }
    //fprintf(TimeTable::report, "delete current: %x\n",TimeTable::current);
    delete TimeTable::current;
    TimeTable::current = NULL;

    TimeTable::VariableValue=0;

    if ( TimeTable::report) {
      fprintf( TimeTable::report, "\n\n");
      // porque estaba comentada la linea?
      //  por el too many open files con el debug [clock] norm activado!
      // que ni siquiera funcionaba tcltk, dejaba guardar, leer, . . .
      //fprintf( "delete report: %x\n",TimeTable::report);
      fclose( TimeTable::report);
      TimeTable::report = NULL;
    }
  }
}

//====================================================================== End
