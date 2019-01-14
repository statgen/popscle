/*
 *  Copyright (C) 2010  Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __COMMANDS_H__
#define __COMMANDS_H__

#include <ctype.h>
#include <stddef.h>
#include <vector>
#include <map>
#include <string>

typedef std::map<std::string,void*> StringMap;
typedef std::map<std::string,void*>::iterator StringMapIterator;

class commandList;  // list of commands

// contents for long command list
struct longCommandList
{
  const char * desc;
  int (*func) (int, char**); 
  const char * help;
};

class command  // represent each command - an abstract class
{
 protected:
  std::string description;  // full name of the option
  int (*func)(int, char**);               // pointer to store actual value
  std::string helpstring;   // detailed description of the option

  static int nameCol;       // length of name column
  static int statusCol;     // length of status column
  static int helpCol;       // length of the help column

  // get option by full string
  //virtual int TranslateExtras(const char * value, const char * extras);
  virtual longCommandList* Translate(const char* value);

  //static bool CheckInteger(const char * value); 
  //static bool CheckDouble(const char * value);

  std::string * errors;
  std::string * messages;

 public:

  // constructor 
  command(const char * desc, int (*f)(int, char**), const char * help = NULL);

  // destructor
  virtual ~command() {}

  // Read argn-th argument and assing values
  //virtual int Read(int argc, char ** argv, int argn) = 0; // {}

  // virtual function which prints the name and values
  virtual void Status() = 0;
  virtual void HelpMessage() = 0;

  // modify nameCol
  static void SetNameLen(int len)
  {
    nameCol = len;
  }

  // modify statusCol
  static void SetStatusLen(int len)
  {
    statusCol = len;
  }

  // set error buffer
  void SetErrorBuffer(std::string & buffer)
  {
    errors = &buffer;
  }

  // set error buffer
  void SetMessageBuffer(std::string & buffer)
  {
    messages = &buffer;
  }

  // function printing warning
  //void error(const char * format, ...);
  //void message(const char * format, ...);

  // give full access to paramList class
  friend class commandList;
};

#define BEGIN_LONG_COMMANDS(array)   longCommandList array[] = {	\
    { NULL,  NULL, NULL},

#define LONG_COMMAND_GROUP(label, help)    { label, NULL, help },
#define LONG_COMMAND(label, funcptr, help) { label, funcptr, help},
#define END_LONG_COMMANDS()                { NULL,  NULL, NULL } };

// long command class
class longCommands : public command
{
 public:
  longCommands(const char * desc, longCommandList * list);  // constructor

  virtual void Status();                     // Print the status
  virtual void HelpMessage();                     // Print the status

  //longCommands * SetPrecision(int precision)   // Set precision for output
  //{
  //  this->precision = precision;
  //
  //  return this;
  //}

 protected:
  std::map<std::string, longCommandList*> index;
  //std::map<std::string, longCommandList*> legacyIndex;

  longCommandList * list;
  int group_len;
  int name_len;
  //int precision;

  virtual longCommandList* Translate(const char* value);  
  //virtual void Translate(const char * value);
  //virtual int TranslateExtras(const char * value, const char * extras);

  void Status(longCommandList * ptr, int & line_len, bool & need_a_comma);
  void HelpMessage(longCommandList * ptr);
};

// List of commands
class commandList {
 protected:
  bool help;
  std::vector<command*> pl; // vector of pointers;
  
 public:
  commandList() : help(false) {}
  virtual ~commandList();

  void Add(command * p);

  // Tries to process all command line arguments
  virtual int Read(int argc, char ** argv, int start = 1);

  // Outputs summary of command switches and settings
  virtual void Status();
  virtual void HelpMessage();

  // Keeps track of warnings generated during command processing
  std::string  errors;
  std::string  messages;
};

#endif // commands.h
