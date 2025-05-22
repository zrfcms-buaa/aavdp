#ifndef __AAVDP_MAIN__
#define __AAVDP_MAIN__
#include <unistd.h>
#include "./MODEL/MODEL.h"
#include "./XRD/XRD.h"
#include "./KED/KED.h"
#include "./DED/DED.h"
#include "./RDF/RDF.h"

#define TYPE_INPUT_NUMBER 10
#define PATH_CHAR_NUMBER 100
#define EXT_CHAR_NUMBER 20

extern void print_version();
extern void print_help();
extern bool is_path_accessible(char path[]);
extern double be_double(char str[]);
extern bool is_parameter(char str[]);
extern bool not_in_switch(char option[], char swich[]);
extern void merge_path(char file_path[], char exts[][EXT_CHAR_NUMBER], int num);
extern void split_path(char name[], char ext[], char file_path[]);

void print_version()
{
    printf("AAVDP Version 1.0.0 (2025.05.22)\n"
           "An integrated command-line program for Atomistic Analyzer of Virtual Diffraction Patterns (AAVDP) of artificial atomistic structures.\n"
           "Copyright[c] 2025-2027, Beihang University by Yan Zhang and Ruifeng Zhang.\n"
           "Please send bugs and suggestions to zrfcms@buaa.edu.cn.\n");
}

void print_help()
{
    printf("The syntax format and rules for AAVDP:\n"
           "    AAVDP <--mode> (inputfile) <-parameter> [value]   \n"
           "    AAVDP --xrd     # X-ray diffraction               \n"
           "    AAVDP --ned     # Neutron diffraction             \n"
           "    AAVDP --ked     # Kinematic electron diffraction\n"
           "    AAVDP --ded     # Dynamical electron diffraction  \n"
           "    AAVDP --kkd     # Kinematic Kikuchi diffraction \n"
           "    AAVDP --dkd     # Dynamical Kikuchi diffraction   \n"
           "    AAVDP --rdf     # Radial distribution function    \n"
           "    AAVDP --ssf     # Static structure factor         \n"
           "Check /man/manual.pdf for more information\n");
}

bool is_path_accessible(char path[])
{
    if(access(path, F_OK)==-1){
        printf("[ERROR] Path %s does not exist", path);
        exit(1);
    }else if(access(path, R_OK)==-1){
        printf("[ERROR] Path %s is not readable", path);
        exit(1);
    }
    return true;
}

double be_double(char str[])
{
    char *err;
    double value=strtod(str, &err);
    if(strlen(err)!='\0')
    {
        printf("[ERROR] Input parameter %s should be float or int type", str);
        exit(1);
    }
    return value;
}

bool is_parameter(char str[])
{
    if(0==strncmp(str, "-", 1)){
        return false;
    }
    return true;
}

bool not_in_switch(char option[], const char swich[])
{
    if((!is_parameter(option))&&strstr(option, swich)==nullptr){
        return true;
    }
    return false;
}

void merge_path(char file_path[], char exts[][EXT_CHAR_NUMBER], int num)
{
    for(int i=0;i<num;i++){
        strcat(file_path, exts[i]);
    }
}

void split_path(char name[], char ext[], char file_path[])
{
    char *ch=strrchr(file_path,'.');
    strcpy(name, file_path);
    strcpy(ext, ch);
    if(ch!=nullptr){
        int pos=ch-file_path;
        name[pos]='\0';
    }
}

#endif