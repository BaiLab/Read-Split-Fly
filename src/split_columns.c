#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>


/*
 Version history... 

 Oct 27, 2015 - fixed bug with filename of output that was causing
 filename to be wrong sometimes, which caused the pipeline to also fail. 
*/

#define SUFFIX ".split1stcolumn"

//D2FC08P1:143:D0KHCACXX:6:1101:1563:1953 2:N:0:-R-68-130/2  

/*
 Updated version of this function...
 pre-condition: len is the index of the first \t in readName. part before first \t is called "the id".
 pre-condition: if this is from a paired-end read, the last 2 characters of the id are either /1 or /2.
 post-condition: if from a paired-end read, last 2 characters of id are put at the beginning of the id instead.
 post-condition: the last 3 '-' characters in readName are replaced with \t, so that later stages will view these as separate fields.
 */
void split_field(char *readName, int len) {
  int i,j, offset;
  const int NUM_DASHES = 3; // number of dashes to expect to need replacing with '\t'

  if (!readName || strlen(readName) <= len) {
    fprintf(stderr,"No field supplied to split_field. aborting\n");
    exit(1);
  }
  //actual output
  char t1,t2;
  if (readName[len-2] == '/' && (readName[len-1] == '1' || readName[len-1] == '2')) { //paired
    //before:
    //D2FC08P1:143:D0KHCACXX:6:1101:1563:1953 2:N:0:-R-68-130/2  
    //let's do some horse-trading...
    //put the /1 or /2 at the beginning of the substring...
    t1 = readName[len - 2];
    t2 = readName[len - 1];
    int numDashesFound = 0;
    //shift everything right
    // removed the 'int' here....
    for (i = (len - 3);i >= 0 && numDashesFound < NUM_DASHES; i--) {
      if (readName[i] == '-') {
        numDashesFound++;
      }
      readName[i+2] = readName[i];
    }
    readName[i+1] = t1;
    readName[i+2] = t2;
  }

  // replace '-' with '\t'
  int numDashesFound = 0;
  for (i = (len - 1);i >= 0 && numDashesFound < NUM_DASHES; i--) {
    if (readName[i] == '-') {
      numDashesFound++;
      readName[i] = '\t'; // replace '-' with '\t'
    }
  }
  
  //after, OR, the /2 is irrelevant
  //D2FC08P1:143:D0KHCACXX:6:1101:1563:1953 2:N:0:/2-R-68-130
  // except that now, - have already been replaced with \t so no need for the next 3 lines,
  // which is good because they assume the length of the read half is at most 99, which it might not be in fact.
  //i = len - 5;
  //readName[i] = '\t';
  //readName[i+2] = '\t';
}

void parse(char *filename) {
  char *line = 0;
  char *readName,*tmp;
  size_t n;
  int pos;
  FILE *in, *out;
  char *outname;

  
  if (!(in = fopen(filename,"r"))) {
    fprintf(stderr, "Could not open %s for reading", filename);
    perror("");
    return;
  }
  
  // bug fixed on Oct 27, 2015 where malloc was not allocating room for the NULL byte
  // at the end of the string. this caused the filename to be wrong sometimes, causing
  // the pipeline to fail sometimes.
  outname = (char *)malloc(sizeof(char) * (strlen(filename)+strlen(SUFFIX)+1));
  memset(outname, 0, sizeof(char) * (strlen(filename)+strlen(SUFFIX)+1));
  sprintf(outname, "%s%s",filename,SUFFIX);
  outname[strlen(filename)+strlen(SUFFIX)] = 0;
  if (!(out = fopen(outname, "w"))) {
    fprintf(stderr,"Could not open %s for writing", outname);
    perror("");
    fclose(in);
    return;
  }

  tmp = 0;
  while((!feof(in)) && (getline(&line,&n,in) != -1)) {
    if (tmp) { free(tmp); tmp=0; }
    tmp = strdup(line);
    char *found = strchr(line, '\t');
    pos = (found - line);
    split_field(tmp, pos);
    fputs(tmp,out);
  }

  fclose(in);
  fclose(out);
}


int main(int argc, char *argv[]) {
  int i;

  if (argc < 2) {
    fprintf(stderr, "usage -- %s <file to split> [additonal files...]\n",argv[0]);
    return 1;
  }

  for (i = 1;i < argc;i++)
    parse(argv[i]);

  return 0;
}
