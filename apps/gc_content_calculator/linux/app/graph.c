#include <stdio.h>
#include <stdlib.h>
#include <glib.h>
#include <cairo.h>
#include <math.h>

#include "structures.h"

void graph(GtkWidget *widget,  TObject *text_struct); // Build a diagram
void draw_callback(GtkWidget *widget, cairo_t *cr, gpointer data); // Drawing a diagram
void save_graph(GtkWidget *widget, gint response_id, TObject *text_struct); // Save the diagram image


void draw_callback(GtkWidget *widget, cairo_t *cr, gpointer data)
{   // Drawing a diagram

    TObject *text_struct = (TObject *)data;

    cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);

    guint width, height;
    GdkRGBA color;

    // We get the dimensions of the widget
    width = gtk_widget_get_allocated_width(widget);
    height = gtk_widget_get_allocated_height(widget) - 30;

    // Setting the background color
    gdk_rgba_parse(&color, "white");
    gdk_cairo_set_source_rgba(cr, &color);
    cairo_paint(cr);

    // Normalize the height of the columns to the height of the widget
    double normalize_factor = 180.0;

    GtkTextBuffer *buffer = gtk_text_view_get_buffer(
        GTK_TEXT_VIEW(text_struct->text_field)
    );
    GtkTextIter start, end;
    gtk_text_buffer_get_start_iter(buffer, &start);
    gtk_text_buffer_get_end_iter(buffer, &end);

    char *text = gtk_text_buffer_get_text(
        buffer, &start, &end, FALSE
    );
    char name[4096];
    char file_header[MAX_HEADER_LENGTH];
    int hour, minute, day, month, year;
    long int nucleotide_counts[5];
    float nucleotide_percent[5];
    long int total_char[4];
    float gc_content_percent;

    int result = sscanf(
        text,
        "Filename: %[^\n]\n\n"
        "File header: %[^\n]\n\n"
        "Date of processing:\n%02d:%02d - %02d.%02d.%d\n\n"
        "A - %ld - %f%%\n"
        "U - %ld - %f%%\n"
        "G - %ld - %f%%\n"
        "C - %ld - %f%%\n"
        "T - %ld - %f%%\n"
        "Total number of characters - %ld\n"
        "Number of ATGCU - %ld\n"
        "Other characters - %ld\n"
        "Number of GC content in characters - %ld\n"
        "Number of GC content in percent - %f%%",
        &name,
        &file_header,
        &hour, &minute, &day, &month, &year,
        &nucleotide_counts[0], &nucleotide_percent[0], // &a_count, &a_percent
        &nucleotide_counts[1], &nucleotide_percent[1], // &u_count, &u_percent
        &nucleotide_counts[2], &nucleotide_percent[2], // &g_count, &g_percent
        &nucleotide_counts[3], &nucleotide_percent[3], // &c_count, &c_percent
        &nucleotide_counts[4], &nucleotide_percent[4], // &t_count, &t_percent
        &total_char[0], // &total_characters
        &total_char[1], // &atgcu
        &total_char[2], // &other_characters
        &total_char[3], // &gc_content_characters
        &gc_content_percent
    );
    if (result != 22) {
        g_free(text);
        char error[] = "Error! Check the data.";
        cairo_set_font_size(cr, 15.0);
        cairo_set_source_rgb(cr, 0, 0, 0);
        cairo_move_to(cr, width / 2 - 100, height / 2 + 30);
        cairo_show_text(cr, g_strdup_printf(error));
    }
    else {
        cairo_set_font_size(cr, 13.0);
        double bar_width = 25.0;

        // Colors for each column
        const char *bar_colors[] = {"red", "green", "blue", "yellow", "orange"};

        double height_y = height - 300;
        char title_text_of_diagram[] = "The number of nucleotides as a percentage:";
        cairo_set_source_rgb(cr, 0, 0, 0);
        cairo_move_to(cr, 450, height_y);
        cairo_show_text(cr, g_strdup_printf(title_text_of_diagram));

        char filename[] = "Filename";
        cairo_set_source_rgb(cr, 0, 0, 0);
        cairo_move_to(cr, 30, height_y - 100 - 25);
        cairo_show_text(cr, g_strdup_printf("%s: %s", filename, name));

        char file_header_name[] = "File header";
        cairo_set_source_rgb(cr, 0, 0, 0);
        cairo_move_to(cr, 30, height_y - 100);
        cairo_show_text(cr, g_strdup_printf("%s: %s", file_header_name, file_header));

        char date[] = "Date of processing";
        cairo_set_source_rgb(cr, 0, 0, 0);
        cairo_move_to(cr, 30, height_y - 75);
        cairo_show_text(
            cr, g_strdup_printf(
                "%s: %0.2d:%0.2d - %0.2d.%0.2d.%0.2d",
                date, hour, minute, day, month, year
            )
        );

        char count_text_of_diagram[] = "Quantity: ";
        cairo_set_source_rgb(cr, 0, 0, 0);
        cairo_move_to(cr, 30, height_y);
        cairo_show_text(cr, g_strdup_printf(count_text_of_diagram));

        for (int i = 0; i < NUMBER_OF_NUCLEOTIDES; i++) {
            height_y = height_y + 25;
            cairo_move_to(cr, 30, height_y);
            cairo_set_source_rgb(cr, 0, 0, 0);
            cairo_show_text(
                cr, g_strdup_printf(
                    "%c - %ld", text_struct->array_of_nucleotides[i],
                    nucleotide_counts[i]
                )
            );
        }

        height_y = height_y + 25;
        const char *other_count_text_of_diagram_name[] = {
            "Total number of characters",
            "Number of ATGCU", "Other characters",
            "Number of GC content in characters"
        };
        size_t other_count_text_of_diagram_name_size = (
            sizeof(other_count_text_of_diagram_name) / sizeof(other_count_text_of_diagram_name[0])
        );
        for (size_t i = 0; i < other_count_text_of_diagram_name_size; i++) {
            height_y = height_y + 25;
            cairo_move_to(cr, 30, height_y);
            cairo_set_source_rgb(cr, 0, 0, 0);
            cairo_show_text(
                cr, g_strdup_printf(
                    "%s - %ld", other_count_text_of_diagram_name[i],
                    total_char[i]
                )
            );
        }

        height_y = height_y + 25;
        cairo_move_to(cr, 30, height_y);
        cairo_show_text(
            cr, g_strdup_printf(
                "Number of GC content in percent - %0.2f%%",
                gc_content_percent
            )
        );

        // We draw columns for each nucleotide
        for (int i = 0; i < NUMBER_OF_NUCLEOTIDES; i++) {
            double x = 450 + ((2 * i + 1) * (bar_width));
            double y = (height - height * nucleotide_percent[i] / normalize_factor);

            // Setting the color for the current column
            gdk_rgba_parse(&color, bar_colors[i]);
            gdk_cairo_set_source_rgba(cr, &color);

            // Draw a column
            cairo_rectangle(cr, x - bar_width / 2, y, bar_width, height - y);
            cairo_fill(cr);

            // Adding a caption under the column
            cairo_move_to(cr, x - 5, height + 17);
            cairo_set_source_rgb(cr, 0, 0, 0); // black color for text
            char *text = g_strdup_printf("%c", text_struct->array_of_nucleotides[i]);
            cairo_show_text(cr, text);
            g_free(text);

            // Add the number of nucleotides near the column
            cairo_move_to(cr, x - 20, y - 10);
            cairo_show_text(cr, g_strdup_printf("%0.2f", nucleotide_percent[i]));

        }
        // Drawing a graph
        cairo_stroke(cr);
    }
}

void graph(GtkWidget *widget,  TObject *text_struct)
{   // Build a diagram

    GtkWidget *dialog = gtk_dialog_new_with_buttons(
        "Diagram", GTK_WINDOW(gtk_widget_get_toplevel(widget)),
        GTK_DIALOG_MODAL,
        "_Save", GTK_RESPONSE_APPLY,
        "_Close", GTK_RESPONSE_CLOSE,
        NULL
    );
    gtk_widget_set_size_request(dialog, 800, 500);
    gtk_window_set_resizable(dialog, FALSE);

    GtkWidget *drawing_area = gtk_drawing_area_new();

    gtk_container_add(GTK_CONTAINER(
        gtk_dialog_get_content_area(GTK_DIALOG(dialog))),
        drawing_area
    );
    gtk_widget_set_size_request(drawing_area, 800, 500);

    g_signal_connect(G_OBJECT(drawing_area), "draw",
        G_CALLBACK(draw_callback), text_struct
    );
    g_signal_connect(G_OBJECT(dialog), "response",
        G_CALLBACK(save_graph), text_struct
    );
    gtk_widget_activate(dialog);
    gtk_widget_show_all(dialog);
    gtk_dialog_run(GTK_DIALOG(dialog));
    gtk_widget_destroy(dialog);
}

void save_graph(GtkWidget *widget, gint response_id, TObject *text_struct)
{   // Save the diagram image

    if (response_id == GTK_RESPONSE_APPLY)
    {
        GtkWidget *dialog = gtk_file_chooser_dialog_new(
            "Save File", NULL,
            GTK_FILE_CHOOSER_ACTION_SAVE,
            "_Cancel", GTK_RESPONSE_CANCEL,
            "_Save", GTK_RESPONSE_ACCEPT,
            NULL
        );

        gtk_file_chooser_set_do_overwrite_confirmation(
            GTK_FILE_CHOOSER(dialog), TRUE
        );

        if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT)
        {
            char *filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));

            // Create an image and save it to a file
            cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, 800, 550);
            cairo_t *cr = cairo_create(surface);
            draw_callback(widget, cr, text_struct); // Calling the surface rendering function
            cairo_surface_write_to_png(surface, filename); // Saving the image to a file
            cairo_destroy(cr);
            cairo_surface_destroy(surface);
            g_free(filename);
        }
        gtk_widget_destroy(dialog);
    }
}
