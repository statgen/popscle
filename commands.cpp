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

#include "commands.h"
#include "Constant.h"
#include "Error.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctype.h>
#include <stdarg.h>

int command::nameCol = 30;
int command::statusCol = 15;
int command::helpCol = 20;

#define CONSOLE_WIDTH 80

// constructor -- assign initial values
command::command(const char * desc, int (*f)(int, char**), const char* help) :
  description(desc == NULL ? "" : desc),
  func(f),
  helpstring(help == NULL ? "" : help)
{
}

// ???
longCommandList* command::Translate(const char *)
{
  return NULL;
}

// constructor of long command
longCommands::longCommands(const char * desc, longCommandList * lst)
  : command(desc, NULL, NULL)
{
  list = lst;

  index.clear();
  group_len = 0;
  name_len = 0;

  longCommandList * ptr = list + 1;  // first command

  while (ptr->desc != NULL) {  // unless it is non-informative 
    if (ptr->func != NULL) {  // if it is a real command
      index[ptr->desc] = ptr;  // map string to param
      int tmp = strlen(ptr->desc);   // if it is just a group
      if (tmp > name_len) name_len = tmp; // modify group_len      
    }
    else {
      int tmp = strlen(ptr->desc);          // if it is just a group
      if (tmp > group_len) group_len = tmp; // modify group_len
    }
    ptr++;
  }

  //precision = 2;  // ???
}

// cstr is [option_name] and extra is value

longCommandList* longCommands::Translate(const char * cstr)
{
  std::map<std::string, longCommandList*>::iterator it = index.find(cstr);

  if ( it != index.end() ) { // keyword was found
    longCommandList* ptr = it->second;
    return ptr;
    //return -1; // ignore?
  }
  else {
    return NULL;
    //return -1; // ignore?
  }
}

// Print the status of command List
void longCommands::Status(longCommandList * ptr, int & line_len, bool & need_a_comma)
{
  std::string state;
  int line_start = group_len ? group_len + 5 : 0;
  //int i;

  if (ptr->func == NULL) {  // if group command, end previous group (if exists) 
    fprintf(stderr, "%s %*s :", need_a_comma ? "\n" : "", group_len + 2, ptr->desc);
    need_a_comma = false;
    line_len = line_start;
  }
  else {                     // otherwise, print argument name and 
    int item_len = 3 + strlen(ptr->desc) + need_a_comma; // + state.size();
    
    if (item_len + line_len > CONSOLE_WIDTH - 2 && line_len > line_start) {
      line_len = line_start;
      fprintf(stderr, "%s\n%*s", need_a_comma ? "," : "", line_len,  "");
      need_a_comma = 0;
      item_len -= 1;
    }

    fprintf(stderr, "%s %s", need_a_comma ? "," : (need_a_comma = true, ""),
	    ptr->desc);

    need_a_comma = true;
    line_len += item_len;
  }
}

// Print the status of command List
void longCommands::HelpMessage(longCommandList * ptr)
{
  //int i;

  if (ptr->func == NULL) {  // if group command, end previous group (if exists) 
    fprintf(stderr, "\n== %s%s%s ==\n", ptr->desc, ptr->help ? " - " : "", ptr->help ? ptr->help : "");
  }
  else {                     // otherwise, print argument name and 
    fprintf(stderr, "   - %-*s%s%s\n", name_len, ptr->desc, " : ", ptr->help ? ptr->help : "");
  }
}

// print the status of the command
void longCommands::HelpMessage()
{
  if (!description.empty() && description[0] != 0)  // group option
    fprintf(stderr, "\n%s:\n", description.c_str());    
    //fprintf(stderr, "\n%s - %s\n", description.c_str(), helpstring.c_str());

  // for the rest of the group, print commands
  for (longCommandList * ptr = list + 1; ptr->desc != NULL; ptr++)  
    HelpMessage(ptr);

  fprintf(stderr, "\n");
}

// print the status of the command
void longCommands::Status()
{
  if (!description.empty() && description[0] != 0)  { // group option
    fprintf(stderr, "\n%s:\n\n", description.c_str());

    fprintf(stderr, "The following commands are available. Ones with \"[]\" are in effect:\n");
  }

  bool need_a_comma = false;
  int  line_len = 0;

  // for the rest of the group, print commands
  for (longCommandList * ptr = list + 1; ptr->desc != NULL; ptr++)  
    Status(ptr, line_len, need_a_comma);

  fprintf(stderr, "\n");
}

// Add command

void commandList::Add(command * p)
{
  p->SetErrorBuffer(errors);
  p->SetMessageBuffer(messages);
  pl.push_back(p);
}

// Read commands from argument
int commandList::Read(int argc, char ** argv, int start)
{
  std::string cmd(argv[start]);

  //notice("pl.size() = %d", pl.size());
  
  for(int32_t i=0; i < (int32_t)pl.size(); ++i) {
    longCommandList* ptr = pl[0]->Translate(argv[start]);
    if ( ptr != NULL ) {
      //notice("%s %x %s", cmd.c_str(), ptr->func, ptr->desc);            
      if ( ( cmd == ptr->desc ) && ( ptr->func != NULL ) ) {
	fprintf(stderr, "[cramore %s] -- %s\n\n", ptr->desc, ptr->help);
	fprintf(stderr, " Copyright (c) 2009-2017 by Hyun Min Kang and Adrian Tan\n");
	fprintf(stderr, " Licensed under the Apache License v2.0 http://www.apache.org/licenses/\n");
	return ptr->func(argc-start, argv+start);
      }
    }
  }
  HelpMessage();
  error("Cannot recognize the command %s. Run %s --help for detailed instruction", cmd.c_str(), argv[0]);
  return 1;  

  /*
  if ( index.find(cmd) == index.end() ) {
    cl.HelpMessage();
    error("Cannot recognize the command %s. Run %s --help for detailed instruction", cmd.c_str(), argv[0]);
    return 1;
  }
  else {
    //argv[start] = argv[0];
    return index[cmd].func(argc-start, argv+start);
    } */
}  

void commandList::HelpMessage()
{
  fprintf(stderr, "\nDetailed instructions of commands are available:\n");

  for (int i=0; i<(int)pl.size(); i++)
    pl[i]->HelpMessage();

  fprintf(stderr, "\n");

  if (errors.size()) {
    error("[E:%s:%d %s] Problems encountered parsing command line:\n\n%s",__FILE__,__LINE__,__FUNCTION__,
	  errors.c_str());
    errors.clear();
  }

  if (messages.size()) {
    ::printf("NOTES:\n%s\n", messages.c_str());
    //messages.clear();
  }
}

// Print the total command list
void commandList::Status()
{
  fprintf(stderr, "\nThe following commands are available:\n");

  for (int i=0; i<(int)pl.size(); i++)
    pl[i]->Status();

  fprintf(stderr, "\nRun with --help for more detailed help messages of each argument.\n");  
  fprintf(stderr, "\n");

  if (errors.size()) {
    error("[E:%s:%d %s] Problems encountered parsing command line:\n\n%s",__FILE__,__LINE__,__FUNCTION__,
	  errors.c_str());
    errors.clear();
  }

  if (messages.size()) {
    ::printf("NOTES:\n%s\n", messages.c_str());
    messages.clear();
  }
}

//

commandList::~commandList()
{
  for (int i = 0; i < (int)pl.size(); i++)
    delete pl[i];
};
