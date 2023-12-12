#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glib.h>
#include <gtk/gtk.h>

#include "structures.h"
#include "util.h"
#include "sqlite3.h"

void last_ten_iter(GtkWidget *widget,  TObject *text_struct); // Output the 10 last iterations
gboolean is_valid_utf8(const gchar *text); // Checking the validity of UTF-8 strings
gchar *convert_to_utf8(const gchar *text, const gchar *from_encoding); // Converting a string to UTF-8
static int callback(void *data, int argc, char **argv, char **azColName); // Callback function for output the last iteration
gboolean update_ui(gpointer user_data); // For output the last iteration
void last_iter(GtkWidget *widget,  TObject *text_struct); // Output the last iteration
void download_db(GtkWidget *widget,  TObject *text_struct); // Download to DB
gboolean delete_db_callback(gpointer user_data); // For delete DB
void delete_db(GtkWidget *widget, TObject *text_struct); // Delete DB
gboolean create_db_callback(gpointer user_data); // For create DB
void create_db(GtkWidget *widget, TObject *text_struct); // Create DB


void last_ten_iter(GtkWidget *widget,  TObject *text_struct)
{   // Output the 10 last iterations

    sqlite3 *db;
    if (sqlite3_open(text_struct->path_to_db, &db) != SQLITE_OK) {
        char db_error[256];
        snprintf(
            db_error, sizeof(db_error),
            "Error: \n%s",
            sqlite3_errmsg(db)
        );
        sqlite3_close(db);
        return;
    }

    char sql[512];
    snprintf(
        sql, sizeof(sql),
        "SELECT a_count, "
        "g_count, "
        "c_count, "
        "t_count, "
        "u_count, "
        "a_percent, "
        "g_percent, "
        "c_percent, "
        "t_percent, "
        "u_percent, "
        "Total_number_of_characters, "
        "Number_of_ATGCU, "
        "Other_characters, "
        "Number_of_GC_content_in_characters, "
        "Number_of_GC_content_in_percent "
        "FROM atgcu "
        "ORDER BY id DESC "
        "LIMIT 10"
    );
}


gboolean is_valid_utf8(const gchar *text)
{   // Checking the validity of UTF-8 strings

    return g_utf8_validate(text, -1, NULL);
}

gchar *convert_to_utf8(const gchar *text, const gchar *from_encoding)
{   // Converting a string to UTF-8

    return g_convert(text, -1, "UTF-8", from_encoding, NULL, NULL, NULL);
}

static int callback(void *data, int argc, char **argv, char **azColName)
{   // Callback function for output the last iteration

    TObject *text_struct = (TObject *)data;

    g_string_truncate(text_struct->result, 0);
    g_string_append_printf(
        text_struct->result,
        "Date of processing:\n%s:%s - %s.%s.%s\n\n"
        "A - %s - %s%%\n"
        "U - %s - %s%%\n"
        "G - %s - %s%%\n"
        "C - %s - %s%%\n"
        "T - %s - %s%%\n"
        "Total number of characters - %s\n"
        "Number of ATGCU - %s\n"
        "Other characters - %s\n"
        "Number of GC content in characters - %s\n"
        "Number of GC content in percent - %s%%\n",
        argv[0], argv[1], argv[2], argv[3], argv[4],
        argv[5], argv[10],
        argv[9], argv[14],
        argv[6], argv[11],
        argv[7], argv[12],
        argv[8], argv[13],
        argv[15],
        argv[16],
        argv[17],
        argv[18],
        argv[19]
    );

    gchar *utf8_count_out = convert_to_utf8(text_struct->result->str, "UTF-8");

    if (!is_valid_utf8(utf8_count_out)) {
        g_free(utf8_count_out);
        return;
    }

    g_string_assign(text_struct->result, utf8_count_out);
    g_free(utf8_count_out);
    g_idle_add(update_ui, text_struct);

    return 0;
}

gboolean update_ui(gpointer user_data)
{   // For output the last iteration

    TObject *text_struct = (TObject *)user_data;
    GtkTextBuffer *buffer = gtk_text_view_get_buffer(
            GTK_TEXT_VIEW(text_struct->text_field)
    );

    if (text_struct->result) {
        gtk_text_buffer_set_text(buffer, text_struct->result->str, -1);
        g_string_free(text_struct->result, TRUE);
        text_struct->result = NULL;
    }
    else {
        char update_error[] = "Error when displaying the result.";
        error_message(&update_error);
    }

    return G_SOURCE_REMOVE;
}


void last_iter(GtkWidget *widget,  TObject *text_struct)
{   // Output the last iteration

    sqlite3 *db;
    if (sqlite3_open(text_struct->path_to_db, &db) != SQLITE_OK) {
        char db_error[256];
        snprintf(
            db_error, sizeof(db_error),
            "Error: \n%s",
            sqlite3_errmsg(db)
        );
        sqlite3_close(db);
        return;
    }

    char sql[512];
    snprintf(
        sql, sizeof(sql),
        "SELECT "
        "hour, "
        "minute, "
        "day, "
        "month,"
        "year, "
        "a_count, "
        "g_count, "
        "c_count, "
        "t_count, "
        "u_count, "
        "a_percent, "
        "g_percent, "
        "c_percent, "
        "t_percent, "
        "u_percent, "
        "Total_number_of_characters, "
        "Number_of_ATGCU, "
        "Other_characters, "
        "Number_of_GC_content_in_characters, "
        "Number_of_GC_content_in_percent "
        "FROM atgcu "
        "ORDER BY id DESC "
        "LIMIT 1"
    );

    text_struct->result = g_string_new(NULL);
    if (sqlite3_exec(db, sql, callback, text_struct, NULL) != SQLITE_OK) {
        char db_error[256];
        snprintf(
            db_error, sizeof(db_error),
            "SQL error: %s\n",
            sqlite3_errmsg(db)
        );
        error_message(&db_error);
        return;
    }
    sqlite3_close(db);
}

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
    int hour, minute, day, month, year;
    long int a_count, g_count, c_count, t_count, u_count;
    float a_percent, g_percent, c_percent, t_percent, u_percent;
    long int total_characters, atgcu, other_characters, gc_content_characters;
    float gc_content_percent;

    int result = sscanf(
        text,
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
        &hour, &minute, &day, &month, &year,
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
    if (result != 20) {
        g_free(text);
        char count_error[256];
        snprintf(
            count_error, sizeof(count_error),
            "Error writing to the table.\n"
            "%d values.\n"
            "Check what you are trying to "
            "write to the database.",
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
            "Error: \n%s",
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
            "hour, minute, day, month, year, "
            "a_count, u_count, g_count, c_count, t_count, "
            "a_percent, u_percent, g_percent, c_percent, t_percent, "
            "Total_number_of_characters, Number_of_ATGCU, Other_characters, "
            "Number_of_GC_content_in_characters, Number_of_GC_content_in_percent"
        ") "
        "VALUES ("
            "%02d, %02d, %02d, %02d, %d, "
            "%d, %d, %d, %d, %d, '%0.2f', '%0.2f', "
            "'%0.2f', '%0.2f', '%0.2f', %d, %d, %d, %d, '%0.2f'"
        ")",
        hour, minute, day, month, year,
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
        return;
    }
    else {
        GtkTextBuffer *buffer = gtk_text_view_get_buffer(
            GTK_TEXT_VIEW(text_struct->text_field)
        );
        gtk_text_buffer_set_text(buffer, "", -1);
        gtk_text_buffer_insert_at_cursor(
            buffer, "OK", -1
        );
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
            "hour TEXT,"
            "minute TEXT,"
            "day TEXT,"
            "month TEXT,"
            "year TEXT,"
            "a_count TEXT,"
            "g_count TEXT,"
            "c_count TEXT,"
            "t_count TEXT,"
            "u_count TEXT,"
            "a_percent TEXT,"
            "g_percent TEXT,"
            "c_percent TEXT,"
            "t_percent TEXT,"
            "u_percent TEXT,"
            "Total_number_of_characters TEXT,"
            "Number_of_ATGCU TEXT,"
            "Other_characters TEXT,"
            "Number_of_GC_content_in_characters TEXT,"
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
