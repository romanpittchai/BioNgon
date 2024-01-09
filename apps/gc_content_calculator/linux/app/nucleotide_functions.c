#include <stdio.h>
#include <stdlib.h>
#include <gtk/gtk.h>
#include <time.h>

#include "structures.h"
#include "util.h"


void read_file(GtkWidget *widget, TObject *text_struct); // Open selected and read file
void gc_percent(TObject *text_struct); // Search GC percent
void percent_nucl(TObject *text_struct); // Search percent of nucleotides
void count_dna(TObject *text_struct); // Counting dna


void read_file(GtkWidget *widget, TObject *text_struct)
{   // Open selected and read file

    if (text_struct->filename != NULL) {
        text_struct->file = fopen(text_struct->filename, "r");
        count_dna(text_struct);

        time_t t;
        struct tm *now;
        time(&t);
        now = localtime(&t);

        percent_nucl(text_struct);
        GtkTextBuffer *buffer = gtk_text_view_get_buffer(
            GTK_TEXT_VIEW(text_struct->text_field)
        );
        gtk_text_buffer_set_text(buffer, "", -1);

        gtk_text_buffer_insert_at_cursor(
            buffer, g_strdup_printf(
                "Filename: %s\n\n",
                text_struct->name_file
            ),
        -1);

        gtk_text_buffer_insert_at_cursor(
            buffer, g_strdup_printf(
                "Date of processing:\n%02d:%02d - %02d.%02d.%d\n\n",
                now->tm_hour, now->tm_min,
                now->tm_mday, now->tm_mon + 1, now->tm_year + 1900
            ),
        -1);

        for (int i = 0; i < NUMBER_OF_NUCLEOTIDES; i++) {
            gtk_text_buffer_insert_at_cursor(
                buffer, g_strdup_printf(
                    "%c - %ld - %0.2f%\n",
                    (char)text_struct->array_of_nucleotides[i],
                    text_struct->array_number_of_nucleotides[i],
                    text_struct->array_percent_of_nucleotides[i]
                ),
            -1);
        }

        gc_percent(text_struct);


        gtk_text_buffer_insert_at_cursor(
            buffer, g_strdup_printf(
                "Total number of characters - %ld\n",
                text_struct->length_dna
            ), -1);
        gtk_text_buffer_insert_at_cursor(
            buffer, g_strdup_printf(
                "Number of ATGCU - %ld\n",
                text_struct->length_ATGCU
            ),
        -1);
        gtk_text_buffer_insert_at_cursor(
            buffer, g_strdup_printf(
                "Other characters - %ld\n",
                text_struct->length_dna
                - text_struct->length_ATGCU
            ),
        -1);
        gtk_text_buffer_insert_at_cursor(
            buffer, g_strdup_printf(
                "Number of GC content in characters - %ld\n",
                text_struct->gc_content
            ),
        -1);
        gtk_text_buffer_insert_at_cursor(
            buffer, g_strdup_printf(
                "Number of GC content in percent - %0.2f%\n",
                text_struct->gc_content_percent
            ),
        -1);

        fclose(text_struct->file);
        text_struct->file = NULL;
        free(text_struct->name_file);
        text_struct->name_file = NULL;
        free(text_struct->filename);
        text_struct->filename = NULL;
        init_struct(text_struct);
    }
    else if (text_struct->filename == NULL) {
        char error_read_file[] = "No file is selected.";
        error_message(&error_read_file, text_struct);
    }
}

void gc_percent(TObject *text_struct)
{   // Search GC percent

    text_struct->gc_content =(
        ((float)text_struct->array_number_of_nucleotides[2]
        + (float)text_struct->array_number_of_nucleotides[3])
    );
    text_struct->gc_content_percent = (
        (((float)text_struct->array_number_of_nucleotides[2]
        + (float)text_struct->array_number_of_nucleotides[3])
        / (float)text_struct->length_dna) * 100
    );
}


void percent_nucl(TObject *text_struct)
{   // Search percent of nucleotides

    for (int i = 0; i < NUMBER_OF_NUCLEOTIDES; i++) {
        text_struct->array_percent_of_nucleotides[i] = (
            ((float)text_struct->array_number_of_nucleotides[i]
            / (float)text_struct->length_dna) * 100
        );
    }
}


void count_dna(TObject *text_struct)
{   // Counting dna

    while ((text_struct->bytesRead = fread(
        text_struct->buffer, 1, sizeof(text_struct->buffer),
        text_struct->file)) > 0) {
        for (size_t i = 0; i < text_struct->bytesRead; i++) {
            if (text_struct->buffer[i] > 64 && text_struct->buffer[i] < 91) {
                (text_struct->length_dna)++;
                for (int j = 0; j < NUMBER_OF_NUCLEOTIDES; j++) {
                    if (text_struct->array_of_nucleotides[j] == text_struct->buffer[i]) {
                        (text_struct->length_ATGCU)++;
                        (text_struct->array_number_of_nucleotides[j])++;
                        break;
                    }
            }   }
        }
    }
}






