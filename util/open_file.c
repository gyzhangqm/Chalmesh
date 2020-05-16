#include "xplot.h"
#include "stupid_compiler.h"

static void wrong_file_type( char info16 ){
  printf("Error: wrong type of file!\n");
  switch( info16 ){

  case '1':
    printf("This is an ascii command file\n");
    break;

  case '2':
    printf("This is a binary curve + mapping file\n");
    break;

  case '3':
    printf("This is a binary overlapping grid file\n");
    break;

  case '4':
    printf("This is an ascii node-point file\n");
    break;

  case '5':
    printf("This is an ascii overlapping grid file\n");
    break;

  default:
    ;

  }
}


FILE *open_binary_file(input_output *io_ptr, char *prompt, char *deflt, 
		     char read_write, int file_type){
  char *token, info[14], template[14], long_name[250], *chalmesh_home;
  FILE *fp;

  token = get_word( io_ptr, prompt, deflt, 1);

  sprintf(template, "#file-type-%i\n", file_type);

/* if the env variable CHALMESH_HOME is not set, look in the current directory */
  if (!(chalmesh_home = getenv("CHALMESH_HOME")))
    sprintf(long_name, "./%s", token);
  else
    sprintf(long_name, "%s/data/%s", chalmesh_home, token);

  if (read_write == 'w') {
    if ((fp = fopen( token, "wb" )) == NULL &&
	(fp = fopen( long_name, "wb" )) == NULL)
      printf("Warning: Failed to open the file `%s'\n", token);
    else if (file_type != 0){
/* write an identifier line at the top line */
      write( fileno(fp), template, 13*sizeof(char) );
    }
  }
  else if (read_write == 'r') {
    if ((fp = fopen( token, "rb" )) == NULL &&
	(fp = fopen( long_name, "rb" )) == NULL)
      printf("Warning: Failed to open the file `%s'\n", token);
    else if (file_type != 0){
      read( fileno(fp), info, 13*sizeof(char) );
      info[13]='\0';
/* error message */
      if (strcmp(template, info)!=0){
	wrong_file_type( info[11] );
	fclose( fp );
	fp = NULL;
      }
    }
  }
  else{
    printf("Unknown read-write status: %c\n", read_write);
    fp = NULL;
  }

  return fp;
}



FILE *open_ascii_file( input_output *io_ptr, char *prompt, char *deflt, 
		      char **file_name, char read_write, int file_type, 
		      int save_on_copy){
  char info[14], template[14], long_name[250], *chalmesh_home;
  FILE *fp=NULL;
  
  *file_name = get_word( io_ptr, prompt, deflt, save_on_copy );
  
  sprintf(template, "#file-type-%i\n", file_type);

/* if the env variable CHALMESH_HOME is not set, look in the current directory */
  if (!(chalmesh_home = getenv("CHALMESH_HOME")))
    sprintf(long_name, "./%s", *file_name);
  else
    sprintf(long_name, "%s/data/%s", chalmesh_home, *file_name);

  if (read_write == 'w'){
    if ((fp = fopen( *file_name, "w" )) == NULL &&
	(fp = fopen( long_name, "w" )) == NULL)
      printf("Warning: Failed to open the file `%s'\n", *file_name);
    else if (file_type != 0){
/* write an identifier line at the top line */
      write( fileno(fp), template, 13*sizeof(char) );
    }
  }
  else if (read_write == 'r'){
    if ((fp = fopen( *file_name, "r" )) == NULL &&
	(fp = fopen( long_name, "r" )) == NULL)
      printf("Warning: Failed to open the file `%s'\n", *file_name);
    else if (file_type != 0){
      read( fileno(fp), info, 13*sizeof(char) );
      info[13]='\0';
      if (strcmp(template, info)!=0){
/* error message */
	wrong_file_type( info[11] );
	fclose( fp );
	fp = NULL;
      }
    }
  }
  return fp;
}


FILE *
open_this_ascii_file( char *file_name, char read_write, int file_type, int quiet){
  char info[14], template[14], *chalmesh_home, long_name[250];
  FILE *fp=NULL;
  
  sprintf(template, "#file-type-%i\n", file_type);

/* if the env variable CHALMESH_HOME is not set, look in the current directory */
  if (!(chalmesh_home = getenv("CHALMESH_HOME")))
    sprintf(long_name, "./%s", file_name);
  else
    sprintf(long_name, "%s/data/%s", chalmesh_home, file_name);

  if (read_write == 'w'){
    if ((fp = fopen( file_name, "w" )) == NULL &&
        (fp = fopen( long_name, "w" )) == NULL){
      if (!quiet) printf("Warning: Failed to open the file `%s'\n", file_name);
    }
    else if (file_type != 0){
/* write an identifier line at the top line */
      write( fileno(fp), template, 13*sizeof(char) );
    }
  }
  else if (read_write == 'r'){
    if ((fp = fopen( file_name, "r" )) == NULL &&
        (fp = fopen( long_name, "r" )) == NULL){
      if (!quiet) printf("Warning: Failed to open the file `%s'\n", file_name);
    }
    else if (file_type != 0){
      read( fileno(fp), info, 13*sizeof(char) );
      info[13]='\0';
      if (strcmp(template, info)!=0){
/* error message */
	if (!quiet) wrong_file_type( info[11] );
	fclose( fp );
	fp = NULL;
      }
    }
  }
  return fp;
}


FILE *open_this_binary_file( char *file_name, int file_type, int quiet){
  char info[14], template[14], *chalmesh_home, long_name[250];
  FILE *fp;
  
  sprintf(template, "#file-type-%i\n", file_type);

/* if the env variable CHALMESH_HOME is not set, look in the current directory */
  if (!(chalmesh_home = getenv("CHALMESH_HOME")))
    sprintf(long_name, "./%s", file_name);
  else
    sprintf(long_name, "%s/data/%s", chalmesh_home, file_name);

  if ((fp = fopen( file_name, "rb" )) == NULL &&
      (fp = fopen( long_name, "rb" )) == NULL){
    if (!quiet) printf("Warning: Failed to open the file `%s'\n", file_name);
  }
  else if (file_type != 0){
    read( fileno(fp), info, 13*sizeof(char) );
    info[13]='\0';
    if (strcmp(template, info)!=0){
/* error message */
      if (!quiet) wrong_file_type( info[11] );
      fclose( fp );
      fp = NULL;
    }
  }
  return fp;
}
