#include "structures.h"

#ifndef UTIL_H_INCLUDED
#define UTIL_H_INCLUDED


// Checking for a table
int table_exists(const char *table_name, const char *db_path, TObject *text_struct);

// Are there any rows in the table
int get_total_records(const char *table_name, const char *db_path, TObject *text_struct);

// Basic initialization of the structure
void init_struct(TObject *text_struct);

// Displays an error window
void error_message(char *error, TObject *text_struct);

// Checking for the absence of an interface element
int check_elem_interface(GObject *element, char *element_name, TObject *text_struct);

// Exit the window
void exit_the_window(GtkWidget *window);


#endif // UTIL_H_INCLUDED
