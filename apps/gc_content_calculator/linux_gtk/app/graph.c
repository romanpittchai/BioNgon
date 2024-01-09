#include <stdio.h>
#include <stdlib.h>
#include <glib.h>
#include <cairo.h>
#include <math.h>

#include "structures.h"

void graph(GtkWidget *widget,  TObject *text_struct); // Build a diagram

// Функция для отрисовки графика
void draw_callback(GtkWidget *widget, cairo_t *cr, gpointer data)
{

    TObject *text_struct = (TObject *)data;

    //if (text_struct->redraw_requested == 1) {
        cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
        cairo_set_font_size(cr, 13.0);

        guint width, height;
        GdkRGBA color;

        // Получаем размеры виджета
        width = gtk_widget_get_allocated_width(widget);
        height = gtk_widget_get_allocated_height(widget) - 30;

        // Задаем цвет фона
        gdk_rgba_parse(&color, "white");
        gdk_cairo_set_source_rgba(cr, &color);
        cairo_paint(cr);

        // Нормализуем высоту столбцов к высоте виджета
        double normalize_factor = 100.0;

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
        int hour, minute, day, month, year;
        long int nucleotide_counts[5];
        float nucleotide_percent[5];
        long int total_char[4];
        float gc_content_percent;

        int result = sscanf(
            text,
            "Filename: %s\n\n"
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
        if (result != 21) {
            g_free(text);
            char count_error[256];
            snprintf(
                count_error, sizeof(count_error),
                "Error! "
                "%d values.\n"
                "Check the data.",
                result
            );
            error_message(&count_error, text_struct);
            return EXIT_FAILURE;
        }
        double bar_width = 25.0;
        // Цвета для каждого столбца
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

        char date[] = "Date of processing";
        cairo_set_source_rgb(cr, 0, 0, 0);
        cairo_move_to(cr, 30, height_y - 100);
        cairo_show_text(cr, g_strdup_printf("%s: %0.2d:%0.2d - %0.2d.%0.2d.%0.2d", date, hour, minute, day, month, year));

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

        // Рисуем столбцы для каждого нуклеотида

        for (int i = 0; i < NUMBER_OF_NUCLEOTIDES; i++) {
            double x = 450 + ((2 * i + 1) * (bar_width));
            double y = (height - height * nucleotide_percent[i] / normalize_factor);

            // Задаем цвет для текущего столбца
            gdk_rgba_parse(&color, bar_colors[i]);
            gdk_cairo_set_source_rgba(cr, &color);

            // Рисуем столбец
            cairo_rectangle(cr, x - bar_width / 2, y, bar_width, height - y);
            cairo_fill(cr);

            // Добавляем подпись под столбцом
            cairo_move_to(cr, x - 5, height + 17);
            cairo_set_source_rgb(cr, 0, 0, 0); // черный цвет для текста
            char *text = g_strdup_printf("%c", text_struct->array_of_nucleotides[i]);
            cairo_show_text(cr, text);
            g_free(text);

            // Добавляем количество нуклеотидов возле столбца
            cairo_move_to(cr, x - 20, y - 10); // -10 - отступ для числа
            cairo_show_text(cr, g_strdup_printf("%0.2f%%", nucleotide_percent[i]));

        }

        // Отрисовываем график

        cairo_stroke(cr);


        text_struct->redraw_requested = 0;
        printf("%d", text_struct->redraw_requested);
    //}



}

void graph(GtkWidget *widget,  TObject *text_struct)
{

    GtkWidget *dialog = gtk_dialog_new_with_buttons(
        "Diagram", GTK_WINDOW(gtk_widget_get_toplevel(widget)),
        GTK_DIALOG_MODAL,
        "_Close", GTK_RESPONSE_CLOSE,
        NULL
    );
    gtk_widget_set_size_request(dialog, 800, 500);
    gtk_window_set_resizable (
        dialog,
        FALSE
    );

    GtkWidget *drawing_area = gtk_drawing_area_new();
    gtk_container_add(GTK_CONTAINER(gtk_dialog_get_content_area(GTK_DIALOG(dialog))), drawing_area);
    gtk_widget_set_size_request(drawing_area, 800, 500);

    g_signal_connect(G_OBJECT(drawing_area), "draw", G_CALLBACK(draw_callback), text_struct);
    gtk_widget_activate(dialog);
    gtk_widget_show_all(dialog);
    gtk_dialog_run(GTK_DIALOG(dialog));
    gtk_widget_destroy(dialog);



    while (gtk_events_pending()) {
        gtk_main_iteration();
    }
    text_struct->redraw_requested = 1;
}
