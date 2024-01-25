#include <stdio.h>
#include <stdlib.h>
#include <glib.h>
#include <gtk/gtk.h>

#include "structures.h"

void about_program(GtkWidget *widget,  TObject *text_struct)
{   // About the program

    char about[4096];
    snprintf(
        about, sizeof(about),
        "V. 1.02\n\n"
        "\tThis is a simple program "
        "for counting nucleotides (GC content calculator)."
        "It is a freely distributed software. "
        "It can be used for any purpose and in any way.\n\n"
        "Files:\n"
        "Open - open a file for processing.\n"
        "Save - save the contents of the text field to a text file.\n"
        "Clear text field - clearing the text field.\n"
        "Exit - exit the program.\n\n"
        "DB:\n"
        "Create DB - create a database.\n"
        "Delete DB - delete a database.\n"
        "Last iteration - "
        "output the last record from the database to the text field.\n"
        "Last 10 iterations - "
        "output the last 10 records from the database to the text field.\n"
        "Upload to CSV - save data from the database to a csv file.\n"
        "Delete last iteration - "
        "delete the last record from the database.\n"
        "Delete 10 last iterations - "
        "delete the last 10 records from the database.\n"
        "Clear DB - completely clear the database.\n\n"
        "Other:\n"
        "Graph - Build a diagram based on the received data.\n"
        "About - information about the program.\n\n"
        "Bottom buttons:\n"
        "Run - Make a calculation.\n"
        "Download to DB - Make an entry in the database.\n"
        "Exit - Exit the program.\n"
    );

    GtkTextBuffer *buffer = gtk_text_view_get_buffer(
        GTK_TEXT_VIEW(text_struct->text_field)
    );
    gtk_text_buffer_set_text(buffer, "", -1);
    gtk_text_buffer_insert_at_cursor(
        buffer, about, -1
    );
}
