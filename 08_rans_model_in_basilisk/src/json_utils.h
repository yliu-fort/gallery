#ifndef JSON_UTILS_H
#define JSON_UTILS_H
#include <stdlib.h>
#include <stdio.h>
#include "cJSON.h"
#include <sys/stat.h>

char* read_file(const char *filename)
{
   char *buffer = NULL;
   int string_size, read_size;
   FILE *handler = fopen(filename, "r");

   if (handler)
   {
       // Seek the last byte of the file
       fseek(handler, 0, SEEK_END);
       // Offset from the first to the last byte, or in other words, filesize
       string_size = ftell(handler);
       // go back to the start of the file
       rewind(handler);

       // Allocate a string that can hold it all
       buffer = (char*) malloc(sizeof(char) * (string_size + 1) );

       // Read it all in one operation
       read_size = fread(buffer, sizeof(char), string_size, handler);

       // fread doesn't set it so put a \0 in the last position
       // and buffer is now officially a string
       buffer[string_size] = '\0';

       if (string_size != read_size)
       {
           // Something went wrong, throw away the memory and set
           // the buffer to NULL
           free(buffer);
           buffer = NULL;
       }

       // Always remember to close the file.
       fclose(handler);
    }

    return buffer;
}

struct rsm_JSONFile_t {
  const char *path;
  time_t mtime;
  cJSON *json;
};

int update_configuration_file(struct rsm_JSONFile_t *jf) {
    struct stat file_stat;
    int err = stat(jf->path, &file_stat);
    if (err != 0) {
        perror(" [file_is_modified] stat");
        exit(err);
        //return 0x2;
    }

    char *string = read_file(jf->path);
    if(string == NULL)
      return 0x3;

    cJSON* new_json = cJSON_Parse(string);
    if(new_json == NULL)
        return 0x4;

    cJSON_Delete(jf->json);
    jf->json = new_json;
    jf->mtime = file_stat.st_mtime;

    free(string);
    return 0x1;
}

int update_configuration_file_when_modified(struct rsm_JSONFile_t *jf) {
    struct stat file_stat;
    int err = stat(jf->path, &file_stat);
    if (err != 0) {
        perror(" [file_is_modified] stat");
        exit(err);
        //return 0x2;
    }
    if(file_stat.st_mtime > jf->mtime)
    {
        return update_configuration_file(jf);
    }
    return 0x0;
}

void rsm_print_json(const cJSON *obj){
    char *string1 = cJSON_Print(obj);
    fprintf (stderr, "%s", string1);
    if (string1)
        free(string1);
}

void rsm_print_jsonfile(struct rsm_JSONFile_t *jf){
    char *string1 = cJSON_Print(jf->json);
    fprintf (stderr, "%s", string1);
    if (string1)
        free(string1);
}

void read_double(double *out, const char* entry, const cJSON *const obj)
{
    const cJSON *val = NULL;
    val = cJSON_GetObjectItemCaseSensitive(obj, entry);
    if (cJSON_IsNumber(val))
    {
        *out = val->valuedouble;
    }
}

void read_int(int *out, const char* entry, const cJSON *const obj)
{
    const cJSON *val = NULL;
    val = cJSON_GetObjectItemCaseSensitive(obj, entry);
    if (cJSON_IsNumber(val))
    {
        *out = val->valueint;
    }
}

void read_bool(bool *out, const char* entry, const cJSON *const obj)
{
    const cJSON *val = NULL;
    val = cJSON_GetObjectItemCaseSensitive(obj, entry);
    if (cJSON_IsBool(val))
    {
        *out = val->valueint;
    }
}

void read_char(char *out, const char* entry, const cJSON *const obj)
{
    const cJSON *val = NULL;
    val = cJSON_GetObjectItemCaseSensitive(obj, entry);
    if (cJSON_IsString(val))
    {
        sprintf(out, "%s", val->valuestring);
    }
}

void read_num_array_1d(double** out, const char* entry, const cJSON *const obj)
{
    double *_out = NULL;
    const cJSON *val = NULL;
    const cJSON *vals = NULL;
    vals = cJSON_GetObjectItemCaseSensitive(obj, entry);
    if(vals == NULL || !cJSON_IsArray(vals)){return;}

    int len = cJSON_GetArraySize(vals);
    _out = calloc(len, sizeof(double));
    int it = 0;
    cJSON_ArrayForEach(val, vals)
    {
        //val = cJSON_GetArrayItem(vals, i);
        if (!cJSON_IsNumber(val)){free(_out);return;}
        _out[it++] = val->valuedouble;
    }

    if(_out != NULL) 
    {
        if(*out != NULL) {free(*out);}
        *out = _out;
    }
}
#endif