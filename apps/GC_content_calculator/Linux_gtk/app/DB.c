#include <gtk/gtk.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "structures.h"
#include "util.h"
#include "sqlite3.h"


void download_db(GtkWidget *widget,  TObject *text_struct)
{   // Download to DB

    GtkTextBuffer *buffer = gtk_text_view_get_buffer(
        GTK_TEXT_VIEW(text_struct->text_field)
    );
    GtkTextIter start, end;
    gtk_text_buffer_get_start_iter(buffer, &start);
    gtk_text_buffer_get_end_iter(buffer, &end);

    char *text = gtk_text_buffer_get_text(
        buffer, &start, &end, FALSE
    );
    long int a_count, g_count, c_count, t_count, u_count;
    float a_percent, g_percent, c_percent, t_percent, u_percent;
    long int total_characters, atgcu, other_characters, gc_content_characters;
    float gc_content_percent;

    int result = sscanf(
        text,
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
        &a_count, &a_percent,
        &u_count, &u_percent,
        &g_count, &g_percent,
        &c_count, &c_percent,
        &t_count, &t_percent,
        &total_characters,
        &atgcu,
        &other_characters,
        &gc_content_characters,
        &gc_content_percent
    );

    if (result != 15) {
        g_free(text);
        char count_error[256];
        snprintf(
            count_error, sizeof(count_error),
            "Error writing to the table. "
            "%d values",
            result
        );
        error_message(&count_error);
        return;
    }

    sqlite3 *db;
    if (sqlite3_open(text_struct->path_to_db, &db) != SQLITE_OK) {
        char db_error[256];
        snprintf(
            db_error, sizeof(db_error),
            "Cannot create database: \n%s",
            sqlite3_errmsg(db)
        );
        sqlite3_close(db);
        g_free(text);
        return;
    }

    char sql[512];
    snprintf(
        sql, sizeof(sql),
        "INSERT INTO atgcu ("
            "a_count, u_count, g_count, c_count, t_count, "
            "a_percent, u_percent, g_percent, c_percent, t_percent, "
            "Total_number_of_characters, Number_of_ATGCU, Other_characters, "
            "Number_of_GC_content_in_characters, Number_of_GC_content_in_percent"
        ") "
        "VALUES ("
            "%d, %d, %d, %d, %d, '%0.2f', '%0.2f', "
            "'%0.2f', '%0.2f', '%0.2f', %d, %d, %d, %d, '%0.2f'"
        ")",
        a_count, u_count, g_count, c_count, t_count,
        a_percent, u_percent, g_percent, c_percent, t_percent,
        total_characters,
        atgcu,
        other_characters,
        gc_content_characters,
        gc_content_percent
    );

    if (sqlite3_exec(db, sql, 0, 0, 0) != SQLITE_OK) {
        char db_error[256];
        snprintf(
            db_error, sizeof(db_error),
            "SQL error: %s\n",
            sqlite3_errmsg(db)
        );
        error_message(&db_error);
    }
    sqlite3_close(db);
    g_free(text);
}

gboolean delete_db_callback(gpointer user_data)
{   // Delete DB

    TObject *text_struct_c = (TObject *)user_data;

    if (remove(text_struct_c->path_to_db) == 0) {
        char db_delete_info[] = "The database has been deleted";
        error_message(&db_delete_info);
    }
    else {
        char db_delete_error[] = "The database could not be deleted";
        error_message(&db_delete_error);
    }

    return G_SOURCE_REMOVE;
}


void delete_db(GtkWidget *widget, TObject *text_struct)
{   // Delete DB

    if (access(text_struct->path_to_db, F_OK) != -1) {
        g_idle_add(delete_db_callback, text_struct);
    }
    else {
        char db_delete_error[] = (
            "There is no such database, "
            "there is nothing to delete"
        );
        error_message(&db_delete_error);
    }
}


gboolean create_db_callback(gpointer user_data)
{   // Create DB

    TObject *text_struct_c = (TObject *)user_data;
    sqlite3 *db;
    int rc = sqlite3_open(
        text_struct_c->path_to_db, &db
    );

    if (rc != SQLITE_OK) {
        char db_error[256];
        snprintf(
            db_error, sizeof(db_error),
            "Cannot create database: \n%s",
            sqlite3_errmsg(db)
        );
        error_message(&db_error);
    }
    else {
        char *err_msg = 0;
        char *sql = (
            "CREATE TABLE atgcu ("
            "id INTEGER PRIMARY KEY,"
            "a_count INTEGER,"
            "g_count INTEGER,"
            "c_count INTEGER,"
            "t_count INTEGER,"
            "u_count INTEGER,"
            "a_percent TEXT,"
            "g_percent TEXT,"
            "c_percent TEXT,"
            "t_percent TEXT,"
            "u_percent TEXT,"
            "Total_number_of_characters INTEGER,"
            "Number_of_ATGCU INTEGER,"
            "Other_characters INTEGER,"
            "Number_of_GC_content_in_characters INTEGER,"
            "Number_of_GC_content_in_percent TEXT)"
        );
        rc = sqlite3_exec(db, sql, 0, 0, &err_msg);
        if (rc != SQLITE_OK) {
            char db_error[256];
            snprintf(
                db_error, sizeof(db_error),
                "The table was not created: \n%s",
                err_msg
            );
            error_message(&db_error);
        }
        char create_db_mess[] = "The database has been created";
        error_message(&create_db_mess);
    }

    sqlite3_close(db);

    return G_SOURCE_REMOVE;
}


void create_db(GtkWidget *widget, TObject *text_struct)
{   // Create DB

    char dir_name_db[] = "../DB";

    if (access(dir_name_db, F_OK) == -1) {
        if (mkdir(dir_name_db, 0777) == -1) {
            char dir_error[] = "The directory has not been created";
            error_message(&dir_error);
        }
    }

    if (access(text_struct->path_to_db, F_OK) != -1) {
        char db_error_alr[] = "The database already exists";
        error_message(&db_error_alr);
    }
    else {
        g_idle_add(create_db_callback, text_struct);
    }
}
