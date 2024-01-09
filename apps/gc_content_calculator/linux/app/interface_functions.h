#include <gtk/gtk.h>

#include "structures.h"

#ifndef INTERFACE_GTK_H_INCLUDED
#define INTERFACE_GTK_H_INCLUDED

// To clear the text field
void clear_text_field(GtkWidget *widget, TObject *text_struct);

// Exit the program
void exit_the_program(GtkWidget *main_window);

// Open the data file
void open_file(GtkWidget *widget, TObject **text_struct);

// Save as
void save_file(GtkWidget *widget, TObject *text_struct);

#endif // INTERFACE_GTK_H_INCLUDED
