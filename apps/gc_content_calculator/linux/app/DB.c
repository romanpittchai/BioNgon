#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glib.h>
#include <gtk/gtk.h>

#include "structures.h"
#include "util.h"
#include "sqlite3.h"

void clear_db(GtkWidget *widget, TObject *text_struct); // Clear DB
void delete_ten_last_iteration(GtkWidget *widget, TObject *text_struct); // Delete last 10 iterations
void delete_last_iteration(GtkWidget *widget, TObject *text_struct); // Delete last iteration
int callback_db(void *data, int argc, char **argv, char **azColName); // For uploading the database to csv
void uploading_database(GtkWidget *widget, TObject *text_struct); // Uploading the database to csv
void append_result_to_string(sqlite3_stmt *stmt, GString *result); // Append result to string
void update_text_view_buffer(GtkWidget *text_view, const gchar *text); // Update text view buffer
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

    if (access(text_struct->path_to_db, F_OK) == -1) {
        char db_dowload_error[] = (
            "The database is missing."
        );
        error_message(&db_dowload_error, text_struct);
        return;
    }

    if (table_exists(text_struct->table_name, text_struct->path_to_db, text_struct) < 1) {
        char db_error[256];
            snprintf(
                db_error, sizeof(db_error),
                "The table is missing -'%s'",
                text_struct->table_name
            );
        error_message(&db_error, text_struct);
        return;
    }

    if (get_total_records(text_struct->table_name, text_struct->path_to_db, text_struct) < 1) {
        char db_dowload_error[] = (
            "There is not a single row in the table."
        );
        error_message(&db_dowload_error, text_struct);
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
        error_message(&db_error, text_struct);
        sqlite3_close(db);
        return;
    }

    char sql[512];
    snprintf(
        sql, sizeof(sql),
        "SELECT "
        "filename, hour,"
        "minute, day, "
        "month,year, "
        "a_count, g_count, "
        "c_count, t_count, "
        "u_count, a_percent, "
        "g_percent, c_percent, "
        "t_percent, u_percent, "
        "Total_number_of_characters, "
        "Number_of_ATGCU, "
        "Other_characters, "
        "Number_of_GC_content_in_characters, "
        "Number_of_GC_content_in_percent "
        "FROM %s "
        "ORDER BY id DESC "
        "LIMIT 10",
        text_struct->table_name
    );

    text_struct->result = g_string_new(NULL);

    sqlite3_stmt *stmt;
    int rc = sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);

    if (rc != SQLITE_OK) {
        char db_error[256];
        snprintf(
            db_error, sizeof(db_error),
            "Error preparing statement:\n %s",
            sqlite3_errmsg(db)
        );
        error_message(&db_error, text_struct);
        sqlite3_close(db);
        return;
    }

    while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
        append_result_to_string(stmt, text_struct->result);
    }

    if (rc != SQLITE_DONE) {
        char db_error[256];
        snprintf(
            db_error, sizeof(db_error),
            "Error executing statement: \n%s",
            sqlite3_errmsg(db)
        );
        error_message(&db_error, text_struct);
        return;
    }

    sqlite3_finalize(stmt);
    sqlite3_close(db);

    if (text_struct->result) {
        GtkWidget *text_view = GTK_TEXT_VIEW(text_struct->text_field);
        update_text_view_buffer(text_view, text_struct->result->str);
        g_string_free(text_struct->result, TRUE);
        text_struct->result = NULL;
    } else {
        char update_error[] = "Error when displaying the result.";
        error_message(&update_error, text_struct);
    }
}

void append_result_to_string(sqlite3_stmt *stmt, GString *result)
{   // Append result to string

    g_string_append_printf(
        result,
        "*********************\n\n"
        "Filename: %s\n\n"
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
        "Number of GC content in percent - %s%%\n\n"
        "*********************\n\n",
        sqlite3_column_text(stmt, 0),
        sqlite3_column_text(stmt, 1), sqlite3_column_text(stmt, 2),
        sqlite3_column_text(stmt, 3), sqlite3_column_text(stmt, 4),
        sqlite3_column_text(stmt, 5), sqlite3_column_text(stmt, 6),
        sqlite3_column_text(stmt, 11), sqlite3_column_text(stmt, 10),
        sqlite3_column_text(stmt, 15), sqlite3_column_text(stmt, 7),
        sqlite3_column_text(stmt, 12), sqlite3_column_text(stmt, 8),
        sqlite3_column_text(stmt, 13), sqlite3_column_text(stmt, 9),
        sqlite3_column_text(stmt, 14),
        sqlite3_column_text(stmt, 16),
        sqlite3_column_text(stmt, 17),
        sqlite3_column_text(stmt, 18),
        sqlite3_column_text(stmt, 19),
        sqlite3_column_text(stmt, 20)
    );
}

void update_text_view_buffer(GtkWidget *text_view, const gchar *text)
{   // Update text view buffer

    GtkTextBuffer *buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));
    gtk_text_buffer_set_text(buffer, text, -1);
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
        "*********************\n\n"
        "Filename: %s\n\n"
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
        "Number of GC content in percent - %s%%\n\n"
        "*********************\n\n",
        argv[0],
        argv[1], argv[2], argv[3], argv[4], argv[5],
        argv[6], argv[11],
        argv[10], argv[15],
        argv[7], argv[12],
        argv[8], argv[13],
        argv[9], argv[14],
        argv[16],
        argv[17],
        argv[18],
        argv[19],
        argv[20]
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
        error_message(&update_error, text_struct);
    }

    return G_SOURCE_REMOVE;
}


void last_iter(GtkWidget *widget,  TObject *text_struct)
{   // Output the last iteration

    if (access(text_struct->path_to_db, F_OK) == -1) {
        char db_dowload_error[] = (
            "The database is missing."
        );
        error_message(&db_dowload_error, text_struct);
        return;
    }

    if (table_exists(text_struct->table_name, text_struct->path_to_db, text_struct) < 1) {
        char db_error[256];
            snprintf(
                db_error, sizeof(db_error),
                "The table is missing -'%s'",
                text_struct->table_name
            );
        error_message(&db_error, text_struct);
        return;
    }

    if (get_total_records(text_struct->table_name, text_struct->path_to_db, text_struct) < 1) {
        char db_dowload_error[] = (
            "There is not a single row in the table."
        );
        error_message(&db_dowload_error, text_struct);
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
        error_message(&db_error, text_struct);
        sqlite3_close(db);
        return;
    }

    char sql[512];
    snprintf(
        sql, sizeof(sql),
        "SELECT "
        "filename, hour, "
        "minute, day, "
        "month, year, "
        "a_count, g_count, "
        "c_count, t_count, "
        "u_count, a_percent, "
        "g_percent, c_percent, "
        "t_percent, u_percent, "
        "Total_number_of_characters, "
        "Number_of_ATGCU, "
        "Other_characters, "
        "Number_of_GC_content_in_characters, "
        "Number_of_GC_content_in_percent "
        "FROM %s "
        "ORDER BY id DESC "
        "LIMIT 1",
        text_struct->table_name
    );

    text_struct->result = g_string_new(NULL);

    if (sqlite3_exec(db, sql, callback, text_struct, NULL) != SQLITE_OK) {
        char db_error[256];
        snprintf(
            db_error, sizeof(db_error),
            "SQL error: %s\n",
            sqlite3_errmsg(db)
        );
        error_message(&db_error, text_struct);
        sqlite3_close(db);
        return;
    }
    sqlite3_close(db);
}

void download_db(GtkWidget *widget,  TObject *text_struct)
{   // Download to DB

    if (access(text_struct->path_to_db, F_OK) == -1) {
        char db_dowload_error[] = (
            "The database is missing."
        );
        error_message(&db_dowload_error, text_struct);
        return;
    }

    if (table_exists(text_struct->table_name, text_struct->path_to_db, text_struct) < 1) {
        char db_error[256];
            snprintf(
                db_error, sizeof(db_error),
                "The table is missing -'%s'",
                text_struct->table_name
            );
        error_message(&db_error, text_struct);
        return;
    }

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
    long int a_count, g_count, c_count, t_count, u_count;
    float a_percent, g_percent, c_percent, t_percent, u_percent;
    long int total_characters, atgcu, other_characters, gc_content_characters;
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
    if (result != 21) {
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
        error_message(&count_error, text_struct);
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
        error_message(&db_error, text_struct);
        sqlite3_close(db);
        g_free(text);
        return;
    }

    char sql[512];
    snprintf(
        sql, sizeof(sql),
        "INSERT INTO %s ("
            "filename, "
            "hour, minute, day, month, year, "
            "a_count, u_count, g_count, c_count, t_count, "
            "a_percent, u_percent, g_percent, c_percent, t_percent, "
            "Total_number_of_characters, Number_of_ATGCU, Other_characters, "
            "Number_of_GC_content_in_characters, Number_of_GC_content_in_percent"
        ") "
        "VALUES ("
            "?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?"
        ")",
        text_struct->table_name
    );

    sqlite3_stmt *stmt;
    int rc = sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);
    if (rc != SQLITE_OK) {
        char db_error[256];
        snprintf(
            db_error, sizeof(db_error),
            "Error prepraring stateement:\n %s",
            sqlite3_errmsg(db)
        );
        error_message(&db_error, text_struct);
        sqlite3_close(db);
        g_free(text);
        return;
    }

    char str_hour[3], str_minute[3], str_day[3], str_month[3], str_year[5];
    char str_a_count[256], str_u_count[256], str_g_count[256], str_c_count[256], str_t_count[256];
    char str_a_percent[256], str_u_percent[256], str_g_percent[256], str_c_percent[256], str_t_percent[256];
    char str_total_characters[256], str_atgcu[256], str_other_characters[256], str_gc_content_characters[256];
    char str_gc_content_percent[256];

    snprintf(str_hour, sizeof(str_hour), "%02d", hour);
    snprintf(str_minute, sizeof(str_minute), "%02d", minute);
    snprintf(str_day, sizeof(str_day), "%02d", day);
    snprintf(str_month, sizeof(str_month), "%02d", month);
    snprintf(str_year, sizeof(str_year), "%d", year);
    snprintf(str_a_count, sizeof(str_a_count), "%ld", a_count);
    snprintf(str_u_count, sizeof(str_u_count), "%ld", u_count);
    snprintf(str_g_count, sizeof(str_g_count), "%ld", g_count);
    snprintf(str_c_count, sizeof(str_c_count), "%ld", c_count);
    snprintf(str_t_count, sizeof(str_t_count), "%ld", t_count);
    snprintf(str_a_percent, sizeof(str_a_percent), "%0.2f", a_percent);
    snprintf(str_u_percent, sizeof(str_u_percent), "%0.2f", u_percent);
    snprintf(str_g_percent, sizeof(str_g_percent), "%0.2f", g_percent);
    snprintf(str_c_percent, sizeof(str_c_percent), "%0.2f", c_percent);
    snprintf(str_t_percent, sizeof(str_t_percent), "%0.2f", t_percent);
    snprintf(str_total_characters, sizeof(str_total_characters), "%ld", total_characters);
    snprintf(str_atgcu, sizeof(str_atgcu), "%ld", atgcu);
    snprintf(str_other_characters, sizeof(str_other_characters), "%ld", other_characters);
    snprintf(str_gc_content_characters, sizeof(str_gc_content_characters), "%ld", gc_content_characters);
    snprintf(str_gc_content_percent, sizeof(str_gc_content_percent), "%0.2f", gc_content_percent);

    sqlite3_bind_text(stmt, 1, name, -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 2, str_hour, -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 3, str_minute, -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 4, str_day, -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 5, str_month, -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 6, str_year, -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 7, str_a_count, -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 8, str_u_count, -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 9, str_g_count, -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 10, str_c_count, -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 11, str_t_count, -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 12, str_a_percent, -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 13, str_u_percent, -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 14, str_g_percent, -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 15, str_c_percent, -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 16, str_t_percent, -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 17, str_total_characters, -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 18, str_atgcu, -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 19, str_other_characters, -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 20, str_gc_content_characters, -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 21, str_gc_content_percent, -1, SQLITE_STATIC);

    rc = sqlite3_step(stmt);
    if (rc != SQLITE_DONE) {
        char db_error[256];
        snprintf(
            db_error, sizeof(db_error),
            "Error executing statement: \n%s",
            sqlite3_errmsg(db)
        );
        error_message(&db_error, text_struct);
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

    sqlite3_finalize(stmt);
    sqlite3_close(db);
    g_free(text);
}

gboolean delete_db_callback(gpointer user_data)
{   // Delete DB

    TObject *text_struct_c = (TObject *)user_data;

    if (remove(text_struct_c->path_to_db) == 0) {
        char db_delete_info[] = "The database has been deleted";
        error_message(&db_delete_info, text_struct_c);
    }
    else {
        char db_delete_error[] = "The database could not be deleted";
        error_message(&db_delete_error, text_struct_c);
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
        error_message(&db_delete_error, text_struct);
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
        error_message(&db_error, text_struct_c);
        return;
    }
    else {
        char *err_msg = 0;
        char sql[512];
        snprintf(
            sql, sizeof(sql),
            "CREATE TABLE %s ("
            "id INTEGER PRIMARY KEY, "
            "filename TEXT, "
            "hour TEXT, minute TEXT, "
            "day TEXT, month TEXT, "
            "year TEXT, a_count TEXT, "
            "g_count TEXT, c_count TEXT, "
            "t_count TEXT, u_count TEXT, "
            "a_percent TEXT, g_percent TEXT, "
            "c_percent TEXT, t_percent TEXT, "
            "u_percent TEXT, "
            "Total_number_of_characters TEXT, "
            "Number_of_ATGCU TEXT, "
            "Other_characters TEXT, "
            "Number_of_GC_content_in_characters TEXT, "
            "Number_of_GC_content_in_percent TEXT)",
            text_struct_c->table_name
        );
        rc = sqlite3_exec(db, sql, 0, 0, &err_msg);
        if (rc != SQLITE_OK) {
            char db_error[256];
            snprintf(
                db_error, sizeof(db_error),
                "The table %s was not created: \n%s",
                text_struct_c->table_name,
                err_msg
            );
            remove(text_struct_c->path_to_db);
            error_message(&db_error, text_struct_c);
            return;
        }
        char create_db_mess[] = "The database has been created";
        error_message(&create_db_mess, text_struct_c);
    }

    sqlite3_close(db);

    return G_SOURCE_REMOVE;
}


void create_db(GtkWidget *widget, TObject *text_struct)
{   // Create DB

    if (access(text_struct->path_to_db_dir, F_OK) == -1) {
        if (mkdir(text_struct->path_to_db_dir, 0777) == -1) {
            char dir_error[] = "The directory has not been created";
            error_message(&dir_error, text_struct);
        }
    }

    if (access(text_struct->path_to_db, F_OK) != -1) {
        char db_error_alr[] = "The database already exists";
        error_message(&db_error_alr, text_struct);
    }
    else {
        g_idle_add(create_db_callback, text_struct);
    }
}

int callback_db(void *data, int argc, char **argv, char **azColName)
{   // For uploading the database to csv

    TObject *text_struct = (TObject *)data;

    if (text_struct->headers_print_csv == 0) {
        for (int i = 0; i < argc; i++) {
            fprintf(text_struct->file_write_data, "%s", azColName[i]);
            if (i < argc - 1) {
                fprintf(text_struct->file_write_data, ", ");
            }
        }
        fprintf(text_struct->file_write_data, "\n\n");
        text_struct->headers_print_csv = 1;
    }

    for (int i = 0; i < argc; i++) {
        fprintf(text_struct->file_write_data, "%s", argv[i]);
        if (i < argc - 1) {
            fprintf(text_struct->file_write_data, ", ");
        }
    }
    fprintf(text_struct->file_write_data, "\n");

    return EXIT_SUCCESS;
}

void uploading_database(GtkWidget *widget, TObject *text_struct)
{   // Uploading the database to csv

    text_struct->headers_print_csv = 0;

    if (access(text_struct->path_to_db, F_OK) == -1) {
        char db_dowload_error[] = (
            "The database is missing."
        );
        error_message(&db_dowload_error, text_struct);
        return;
    }

    if (table_exists(text_struct->table_name, text_struct->path_to_db, text_struct) < 1) {
        char db_error[256];
            snprintf(
                db_error, sizeof(db_error),
                "The table is missing -'%s'",
                text_struct->table_name
            );
        error_message(&db_error, text_struct);
        return;
    }

    if (get_total_records(text_struct->table_name, text_struct->path_to_db, text_struct) < 1) {
        char db_error[] = (
            "There is not a single row in the table."
        );
        error_message(&db_error, text_struct);
        return;
    }
    else {
        GtkWidget *dialog;
        GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_SAVE;
        gint res;

        dialog = gtk_file_chooser_dialog_new(
            "Save CSV File", GTK_WINDOW(gtk_widget_get_toplevel(widget)),
            action,
            "Save",
            GTK_RESPONSE_ACCEPT,
            "Cancel",
            GTK_RESPONSE_CANCEL,
            NULL
        );

        res = gtk_dialog_run(GTK_DIALOG(dialog));

        if (res == GTK_RESPONSE_ACCEPT) {
            text_struct->filename_save_as = gtk_file_chooser_get_filename(
                GTK_FILE_CHOOSER(dialog)
            );
        }
        else if (res == GTK_RESPONSE_CANCEL) {
            text_struct->filename_save_as == NULL;
            gtk_widget_destroy(dialog);
            return;
        }

        if (text_struct->filename_save_as != NULL) {
            text_struct->file_write_data = fopen(
                text_struct->filename_save_as, "w+t"
            );
            if (text_struct->file_write_data != NULL) {
                sqlite3 *db;
                if (sqlite3_open(text_struct->path_to_db, &db) != SQLITE_OK) {
                    char db_error[256];
                    snprintf(
                        db_error, sizeof(db_error),
                        "Error: \n%s",
                        sqlite3_errmsg(db)
                    );
                    error_message(&db_error, text_struct);
                    fclose(text_struct->file_write_data);
                    sqlite3_close(db);
                    return;
                }

                char sql[512];
                snprintf(
                    sql, sizeof(sql),
                    "SELECT * FROM %s",
                    text_struct->table_name
                );
                int rc;
                rc = sqlite3_exec(db, sql, callback_db, text_struct, NULL);
                if (rc != SQLITE_OK) {
                    char db_error[256];
                    snprintf(
                        db_error, sizeof(db_error),
                        "Error: \n%s",
                        sqlite3_errmsg(db)
                    );
                    error_message(&db_error, text_struct);
                    fclose(text_struct->file_write_data);
                    sqlite3_close(db);
                    return;
                }
            }
            fclose(text_struct->file_write_data);
        }
        gtk_widget_destroy(dialog);
    }
}

void delete_last_iteration(GtkWidget *widget, TObject *text_struct)
{   // Delete last iteration

    if (access(text_struct->path_to_db, F_OK) == -1) {
        char db_dowload_error[] = (
            "The database is missing."
        );
        error_message(&db_dowload_error, text_struct);
        return;
    }

    if (table_exists(text_struct->table_name, text_struct->path_to_db, text_struct) < 1) {
        char db_error[256];
            snprintf(
                db_error, sizeof(db_error),
                "The table is missing -'%s'",
                text_struct->table_name
            );
        error_message(&db_error, text_struct);
        return;
    }

    if (get_total_records(text_struct->table_name, text_struct->path_to_db, text_struct) < 1) {
        char db_error[] = (
            "There is not a single row in the table."
        );
        error_message(&db_error, text_struct);
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
        error_message(&db_error, text_struct);
        sqlite3_close(db);
        return;
    }

    char sql[512];
    snprintf(
        sql, sizeof(sql),
        "DELETE FROM %s "
        "WHERE id = ("
        "SELECT id FROM "
        "%s ORDER BY id "
        "DESC LIMIT 1)",
        text_struct->table_name,
        text_struct->table_name
    );

    if (sqlite3_exec(db, sql, NULL, NULL, NULL) != SQLITE_OK) {
        char db_error[256];
        snprintf(
            db_error, sizeof(db_error),
            "SQL error: %s\n",
            sqlite3_errmsg(db)
        );
        error_message(&db_error, text_struct);
        sqlite3_close(db);
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
}

void delete_ten_last_iteration(GtkWidget *widget, TObject *text_struct)
{   // Delete last 10 iterations

    if (access(text_struct->path_to_db, F_OK) == -1) {
        char db_dowload_error[] = (
            "The database is missing."
        );
        error_message(&db_dowload_error, text_struct);
        return;
    }

    if (table_exists(text_struct->table_name, text_struct->path_to_db, text_struct) < 1) {
        char db_error[256];
            snprintf(
                db_error, sizeof(db_error),
                "The table is missing -'%s'",
                text_struct->table_name
            );
        error_message(&db_error, text_struct);
        return;
    }

    if (get_total_records(text_struct->table_name, text_struct->path_to_db, text_struct) < 1) {
        char db_error[] = (
            "There is not a single row in the table."
        );
        error_message(&db_error, text_struct);
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
        error_message(&db_error, text_struct);
        sqlite3_close(db);
        return;
    }

    char sql[512];
    snprintf(
        sql, sizeof(sql),
        "DELETE FROM %s "
        "WHERE id IN ("
        "SELECT id FROM "
        "%s ORDER BY id "
        "DESC LIMIT 10)",
        text_struct->table_name,
        text_struct->table_name
    );

    if (sqlite3_exec(db, sql, NULL, NULL, NULL) != SQLITE_OK) {
        char db_error[256];
        snprintf(
            db_error, sizeof(db_error),
            "SQL error: %s\n",
            sqlite3_errmsg(db)
        );
        error_message(&db_error, text_struct);
        sqlite3_close(db);
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
}

void clear_db(GtkWidget *widget, TObject *text_struct)
{   // Clear DB

    if (access(text_struct->path_to_db, F_OK) == -1) {
        char db_dowload_error[] = (
            "The database is missing."
        );
        error_message(&db_dowload_error, text_struct);
        return;
    }

    if (table_exists(text_struct->table_name, text_struct->path_to_db, text_struct) < 1) {
        char db_error[256];
            snprintf(
                db_error, sizeof(db_error),
                "The table is missing -'%s'",
                text_struct->table_name
            );
        error_message(&db_error, text_struct);
        return;
    }

    if (get_total_records(text_struct->table_name, text_struct->path_to_db, text_struct) < 1) {
        char db_error[] = (
            "There is not a single row in the table."
        );
        error_message(&db_error, text_struct);
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
        error_message(&db_error, text_struct);
        sqlite3_close(db);
        return;
    }

    char sql[512];
    snprintf(
        sql, sizeof(sql),
        "DELETE FROM %s",
        text_struct->table_name
    );

    if (sqlite3_exec(db, sql, NULL, NULL, NULL) != SQLITE_OK) {
        char db_error[256];
        snprintf(
            db_error, sizeof(db_error),
            "SQL error: %s\n",
            sqlite3_errmsg(db)
        );
        error_message(&db_error, text_struct);
        sqlite3_close(db);
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
}
