/*  Author: Siva Dharman Naidu <sdharmannaidu@sycamores.indstate.edu> */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<dirent.h>
#include<libgen.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<unistd.h>

#define B 1024000
#define MAX 100
#define TYPE_SPLSEQ ".splSeq"
#define TYPE_RESULTS ".results"
#define TYPE_NOTNOVEL ".splSeqNotNovel"
#define TYPE_NOVEL ".splSeqNovel"
#define FASTA ".FASTA"
#define newFASTA ".new"
#define DUP ".dup"
#define DELM "_"
#define ANG ">"

#define HASHMASK 4194303
#define HASHTABLE_SIZE 4194304


int files_count;
char *filetype;

FILE *fastaFile,*newFASTAFile,*dupFile,*seqFile;

typedef struct node
{
  int hashkey;
  char *f;
  char *s;
  struct node *link;
}
NODE;
NODE *hashtable[HASHTABLE_SIZE];

int endswith(char *string,char *str );
FILE* openfile(char *name, char *mode);
void createfolder(char *name);
char **getfilesinfolder(char *path,char opts);
void readfiles(char *inputfolder,char **names, char *bname);
int hash(char *s);
int readfile(char *f);
void inittable();
void insert(char *f,char *s);
NODE* search(char *k);
void makeFASTASeq(char *path,char *bname);

int main(int args, char *argv[])
{
  char *inputfolder;
  char *foldername;
  char *type_options;
  char **resultfiles;
  pid_t pid;

  if(args == 3)
  {
    foldername = basename(argv[1]);
		//foldername = strdup("blast");
    inputfolder = argv[1];
    type_options = argv[2];
    
    resultfiles = getfilesinfolder(inputfolder,type_options[0]);
    if(files_count > 0)
    {
      printf("1. Found %d result files \n",files_count);
      inittable();
      createfolder(foldername);
			//createfolder("blast");
      readfiles(inputfolder,resultfiles,foldername);
      //readfiles(inputfolder,resultfiles,basename(inputfolder));
      /*pid = fork();   // If we do this, then we delete the database before we run queries......
      if( pid == 0 ) {
        if( !access("/usr/bin/rm", X_OK) ) {
          execl("/usr/bin/rm", "rm", "-fr", foldername, (char *) NULL);
        }
        else {
          execlp("rm", "rm", "-fr", foldername, (char *) NULL);
        }
      }*/
    }
  }
  else
  {
    printf("usage example: makefasta encode_folder_location type_options \n type_options A - All files, B -  results files only, C - .splSeqNotNovel, D - .splSeqNovel\n");
    exit(0);
  }
}

void createfolder(char *name)
{
  struct stat st = {0};
  if(stat(name,&st) == -1)
  {
     mkdir(name, 0700);
     printf("2. Created folder: %s\n",name);
  }
}
FILE* openfile(char *name, char *mode)
{
  FILE *tmp;
  tmp = fopen(name,mode);
  if(tmp == NULL)
  {
    printf("Unable to open the file :%s in %s mode\n",name,mode);
    exit(0);
  }
  return tmp;
}
void closefile(FILE *fd)
{
  fclose(fd);
}
char *replace_s(char *str, char *orig, char *rep)
{
  static char buffer[4096];
  char *p;
  if(!(p = strstr(str,orig)))
     return str;
  strncpy(buffer,str,p-str);
  buffer[p-str] = '\0';
  sprintf(buffer+(p-str),"%s%s",rep,p+strlen(orig));

  return buffer;
}
int endswith(char *string,char *str )
{
  string = strrchr(string, '.');
  if(string != NULL )
    return strcmp(string, str);

  return -1 ;
}
char **getfilesinfolder(char *path,char opts)
{
  char **files;
  struct dirent *dr;
  DIR *d = opendir(path);
  if(d)
  { 
    files = (char **)malloc(MAX * sizeof(char *));
    while((dr = readdir(d))!=NULL)
    {
      char *tmp = strdup(dr->d_name);
      switch(opts)
      {
        case 'A':
             if(endswith(tmp,TYPE_RESULTS)==0)
             {
               files[files_count] = strdup(tmp);
               files_count++;
               filetype = TYPE_RESULTS;
             }
               break;
        case 'B':
             if(endswith(tmp,TYPE_SPLSEQ)==0)
             {
               files[files_count] = strdup(tmp);
               files_count++;
               filetype= TYPE_SPLSEQ;
             }
               break;
        case 'C':
            if(endswith(tmp,TYPE_NOTNOVEL)==0)             
            {
              files[files_count] = strdup(tmp);
              files_count++;
              filetype = TYPE_NOTNOVEL;
            };
               break;
        case 'D':
            if(endswith(tmp,TYPE_NOVEL)==0)
            {
               files[files_count] = strdup(tmp);
               files_count++; 
               filetype = TYPE_NOVEL;
            }
              break;
       default:
            printf("Invalid file type options\n");                         
              break;
      }
   }
 }
 return files;
}

//readfiles(inputfolder,resultfiles,foldername);
void readfiles(char *inputfolder,char **names, char *bname)
{ 
  int i,len,tot = 0;
  char *p;
  char cmd[1024];
  char *tmp = malloc(strlen(bname)+strlen("/")+strlen(filetype)+strlen(bname)+strlen(FASTA)+1);
  len = strlen(inputfolder);
  sprintf(tmp,"%s/%s%s%s",bname,bname,filetype,FASTA);

  char *tmp1 = malloc(strlen(tmp)+strlen(DUP)+1);
  sprintf(tmp1,"%s%s",tmp,DUP);
  char *tmp2 = malloc(strlen(tmp)+strlen(newFASTA)+1);
  sprintf(tmp2,"%s%s",tmp,newFASTA);

  printf("tmp = %s, tmp1 = %s\n",tmp,tmp1);

  fastaFile = openfile(tmp,"w");

  for(i = 0;i<files_count;i++)
  {
    p = malloc(len+strlen("/")+strlen(names[i])+1);
    sprintf(p,"%s/%s",inputfolder,names[i]);
    printf("File: %s\n",p);
    makeFASTASeq((char*)strdup(p),(char*)strdup(bname));        
  }
  
  closefile(fastaFile);
  newFASTAFile = openfile(tmp2,"w");
  dupFile = openfile(tmp1,"w");
  tot =  readfile( (char*)strdup(tmp));
  printf("%s Contains Total No of Lines = %d\n",tmp,tot);
  
  sprintf(cmd,"makeblastdb -in %s -parse_seqids -dbtype nucl > %s.log",tmp2,tmp2);
  system(cmd);
  closefile(dupFile);
  closefile(newFASTAFile);
  /*if(tot > 0)
    makeblastdb();
  else
    fprintf(logFile,"%s does not have any spliced sequences\n",tmp);*/

}
void makeFASTASeq(char *path,char *bname)
{
  int ln = 0;
  int col = 0;
  char b[B];
  char *p,*s1,*s2,*s3,*s4,*s5,*s6;
  char *col0 = malloc(strlen(ANG)+strlen(bname)+strlen(DELM)+1);
  sprintf(col0,"%s%s%s",ANG,bname,DELM);
    seqFile = openfile(path,"r");
  
  while(fgets(b,B,seqFile)!=NULL)
  {
    if(ln == 23)
    {
      b[strlen(b)-1] = '\0';
      p = strtok(b,"\t");
      while(p!=NULL)
      {
				// s1 = >groupname_gene
				// 0 : gene
        if(col == 0)
        {
          s1 = malloc(strlen(col0)+strlen(p)+strlen(DELM)+1);
          sprintf(s1,"%s%s_",col0,p);
        }
				// 1 : chromosome
        else if(col == 1)
        {
          s2 = malloc(strlen(s1)+strlen(p)+strlen(DELM)+1);
          sprintf(s2,"%s%s_",s1,p);
        }
				// 2 : supports
        else if(col == 4)
        {
          s3 = malloc(strlen(s2)+strlen(p)+strlen(DELM)+1);
          sprintf(s3,"%s%s_",s2,p);
        }
				// 3 : splice_length
        else if(col == 5)
        {
          s4 = malloc(strlen(s3)+strlen(p)+strlen(DELM)+1);
          sprintf(s4,"%s%s_",s3,p);
        }
				// 4 : range_supporting_reads
        else if(col == 6)
        {
          s5 = malloc(strlen(s4)+strlen(p)+strlen(DELM)+1);
          p = replace_s(p,"]--[",DELM);
          sprintf(s5,"%s%s",s4,p);
        }
				// Sequence... - Aaron
        else if(col == 9)
        {
          s6 = malloc(strlen(s5)+strlen(p)+strlen("\n")+1);
          sprintf(s6,"%s\n%s",s5,p);
        }
        p = strtok(NULL,"\t");
        col++;
      }
      col = 0;
      fprintf(fastaFile,"%s\n",s6);
     }
    else
      ln++;
  }
  closefile(seqFile);
}
void inittable()
{
  int i;
  for(i = 0;i<HASHTABLE_SIZE;i++)
    hashtable[i] = NULL;
}

NODE* search(char *k)
{
  NODE *entry = hashtable[hash(k)];
  while(entry)
  {
    if(strcmp(entry->f,k) == 0)
       return entry;
    entry = entry->link;
  }
  return 0;
}
void insert(char *f,char *s)
{
  if(search(f)!=0)
    fprintf(dupFile,"%s\n%s\n",f,s);
  else
  {
    int h = hash(f);
    NODE *entry = malloc(sizeof(NODE));
    entry->f = strdup(f);
    entry->s = strdup(s);
    entry->link = hashtable[h];
    hashtable[h] = entry;
    fprintf(newFASTAFile,"%s\n%s\n",f,s);
  }
}

int readfile(char *f)
{
  FILE *fd;
  char b[B];
  char *fline,*sline;
  int i = 0,j = 0;
  int flag = 0;
  int tot = 0;
  fd = openfile(f,"r");

  while(fgets(b,B,fd))
  {
    b[strlen(b)-1] = 0;
    if(b[0] == '>')
    {
      flag = 1;
      fline = strdup((char*)b);
      j++;
    }
    else
    {
      if(flag == 1)
      {
        sline = strdup(b);
        flag= 0;
        j++;
      }
    }
    if(j >= 2)
    {
      /*printf("%s, %s\n",sline,fline);*/
      insert(fline,sline);
      j=0;
    }
    i++;
    tot++;
  }
  closefile(fd);
  return tot;
}
int hash(char *s)
{
  int i;
  int hashval = 0;
  for(i = 0;s[i]!=0;i++)
    hashval = (25 * hashval + s[i]-'a') & HASHMASK;
  return hashval;
}
