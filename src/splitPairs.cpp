/*
  File:        splitPairs.cpp

  Copyright 2016 Jeff Kinne, Yongsheng Bai, Brandon Donham, Aaron Cox, Feng Jiang.
  Permission to use for academic, non-profit purposes is granted, in which case this
  copyright notice should be maintained and original authors acknowledged.

  Author:      Jeff Kinne, jkinne@cs.indstate.edu

  Contents:    Program to take alignments of pieces of unaligned reads and determine
               which could be "split" or "matched" pairs - indicating that the
               unaligned read resulted from a splice.  Also, determine which
               matched pairs support each other (resulted from the same splice junction).

  To compile: g++ splitPairs.cpp -o sp -O4 -std=c++11

  To run:     ./sp options.txt

              Where options.txt is an options file.  If the program is run
              with no command-line arguments it by default processes
              RSW_test.txt with the same parameters as the scripts downloaded
              from Yongsheng Bai's website.  When you run the program it
              prints which files output is written to.  See readOptionsFromFile
              function for the order of the parameters in the options file.

Modification history...  

3/17/2016   - update the algorithm for selecting supporting reads to be more 
             memory-efficient.  This includes writing to the results files as 
      the program runs.  Also, update the .splitPairs format to take
      up less space.

1/6/2016   - bug fixes in a few of the .c files that were not allocating enough space for filenames, causing pipeline to crash sometimes.  Replace a few more %i with %li in .c(pp) files to get rid of compile warnings.
11/24/2015 - Modify RSW.h and this file to print out the actual sequence for split pairs.
             This will appear in the output in the .results files.  
           - Update comments in this file and RSW.h.  Replace %i with %li where required
             for 64 bit systems (when printing results of time(NULL) or vector.size())

8/15/2015 - version 1.0.0 of RSR on github.

* sp4 - 7/30/2015
* 7/30/2015 - fix bug in computing supporting reads that appeared in sp3.
* 7/30/2015 - fixed bug where compare_data was being used in main loop that computes
              matched pairs; compare_data was written to be used in sorting, but was
       being used looking for exact match between lines.  Just stopped using
       it in the main loop and put in the if tests there (to avoid repeating
       this mistake in the future).
* 7/30/2015 - had been working on using openmp to parallelize the main loops.  was
              getting seg fault for unknown reason, so commented out all #pragma
       lines that were part of that effort.
* sp3 - 7/29/2015
* 7/29/2015 - Jeff Kinne - fixed bug that was reporting the number of supporting reads
              as one two small for each junction.
* 7/9/2015 - Jeff Kinne - store strings in stringTable to avoid duplicating the 
             strings.  Should save memory and speed up (because compare string pointers
      rather than strcmp).
* 7/8/2015 - started tracking mod history...
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include <limits.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <time.h>
#include <sys/resource.h>
//#include <omp.h>
using namespace std;

#include "RSW.h"


// parameters input from options file
int maxDistance;  // max difference in aligned pieces to be considered matched pair

  // this field no longer used, software automatically detects read length for each read in an earlier part of the pipeline
//int sampleLength; // length of reads in this data set

int minSpliceLength; // minimum size of splice to consider
int supportPosTolerance; // "buffer boundary" parameter for determining supporting reads
int minSupportingReads;  // only output junctions with at least this many supporting reads
char const *sampleDataFile; // input file name
char const *refFlatFile;    // refFlat file of gene locations
char const *refFlatBoundaryFile; // refFlat file of known intron/extron boundaries
char const *resultsBaseName;     // base file name used for output file names

char buff[MAX_STR_LEN];

char options[MAX_STR_LEN];

time_t beginTime, endTime;  // for keeping track of running time of program

// string table is used to reduce memory usage of program.  for any string
// we need we store it in the string table and then only use the pointer to the
// string any time we need it.  So two records from the input file with the same
// read id won't take up twice as much memory.  See reading of data file for
// how this works.
unordered_set<string> stringTable;

// used in printing of results to avoid printing multiple junction sites
// multiple times
//unordered_set<const char *> readIdsReported;

// data read into the program
vector<struct RSW> data; // from input file of alignments of pieces of unaligned reads
vector<struct RSW_Known> data_known; // from refFlat
vector<struct RSW_Boundaries> data_boundaries; // from refFlat intron/extron boundaries

unordered_map<const char *, RSW_half_data> data_halves; // stores information about max/min length seen from each half of an id
string halfStatsString=""; // computed once all data is read in, then printed later.

const char unfound_string[100] = "UNFOUND_";   // gene name of any junction outside of genes

vector<RSW_splice *> data_splice; // used to store possible jucntions, see RSW.h for RSW_splice definition

int numDifferentReads; // counter...

// function not currently used
int compute_hash(RSW *d) {
  int h = 1, i;
  for(i=0; d->id[i] != '\0'; i++) h = h * d->id[i] % SHRT_MAX;
  for(i=0; d->chromosome[i] != '\0'; i++) h = h * d->chromosome[i] % SHRT_MAX;
  h = h * d->direction % SHRT_MAX;
  return h;
}


char sLine[MAX_LINE+1];
char temp[MAX_LINE+1];

/*
  Function: read_data, read in data file into data vector

  Parameters: filename - file to open and read

  Note: if file is .gz or .lrz then attempt to unzip before reading.  This will
  only work if gunzip and/or lrunzip can be run from the current directory.
*/
void read_data(const char *filename) {
  // open file for reading (from pipe if trying to unzip)
  FILE *f;
  int len = strlen(filename);
  if (len > 3 && strcmp(filename+len-3,".gz")==0) {
    sprintf(temp,"gunzip -c %s", filename);
    f = popen(temp, "r");
  }
  else if (len > 4 && strcmp(filename+len-4,".lrz")==0) {
    sprintf(temp,"cat %s | ./lrunzip", filename);
    f = popen(temp, "r");
  }
  else f = fopen(filename, "r");
  if (f == NULL) {printf("Error reading from file %s\n", filename); exit(0); }

  // read data file one line at a time.
  int result=1;
  while (result > 0) {
    int i; char dir;
    RSW *r = new RSW;
    result = get_line(f, sLine, MAX_LINE);
    if (result < 0) {
      printf("Error reading data file %s, line exceeded %i characters.\n", filename, MAX_LINE);
      delete r;
      break;
    }

    // break line into fields, separated by tab
    char *tempA = strtok(sLine, "\t");
    i=0;
    while (tempA != NULL) {
      string temp = tempA; temp.shrink_to_fit();
      pair<unordered_set<string>::iterator,bool> result;
      switch (i) {
      case 0: // id of read
        // put into string table and store pointer to string.  note
        // that insert just returns a pointer if the string already was in the string table.
        result = stringTable.insert(temp);
        r->id = (result.first)->c_str();
        break;
      case 1: // side
        r->side = temp[0];
        break;
      case 2: // length of piece
        r->length = atoi(temp.c_str());
        break;
      case 3: // total length of read this piece is in
        r->totalReadLength = atoi(temp.c_str());
        break;
      case 4: // direction
        r->direction = temp[0];
        break;
      case 5: // chromosome
        result = stringTable.insert(temp);
        r->chromosome = (result.first)->c_str();
        break;
      case 6: // position
        r->position = atol(temp.c_str());
        break;
      case 9: // count, unused currently
        r->count = atoi(temp.c_str());
        break;
      }
      i++;
      tempA = strtok(NULL,"\t");
      if (tempA == NULL) break;
    }
    if (i < 10) {
      delete r; 
      break;
    }

    r->hash = compute_hash(r); // not currently used

    // update RSW_half information ...

    // look to see if we've seen this id and half before.
    // id is pointer to string with id + which side.
    pair<unordered_set<string>::iterator,bool> res;
    string temp;
    temp = r->id; temp += r->side; temp.shrink_to_fit();
    res = stringTable.insert(temp);
    const char * key = (res.first)->c_str();
    
    temp = r->id; temp += toupper(r->side) == 'L' ? 'R' : 'L';
    res = stringTable.insert(temp);
    const char * otherKey = (res.first)->c_str();
    
    auto h_find = data_halves.find(key);
    
    if (h_find == data_halves.end()) { // if new, insert
      RSW_half_data hd; 
      hd.minLength = hd.maxLength = r->length;
      data_halves.insert({key, hd});
    }
    else { // if not new, update as appropriate
      if (r->length < h_find->second.minLength) {
        h_find->second.minLength = r->length;
      }
      if (r->length > h_find->second.maxLength) {
        h_find->second.maxLength = r->length;
      }
    }
    r->halfKey = key;
    r->otherHalfKey = otherKey;

    // what position on this would be at the split
    if (r->side == 'L' && r->direction == '+' ||
        r->side == 'R' && r->direction == '-') {
      r->splitPos = r->position + r->length;
    }
    else {
      r->splitPos = r->position;
    }
        
    data.push_back(*r);        // save into vector
    delete r;
    //if (data.size() >= 15000000) break; // cut off early, for debugging to prevent program from running for too long.
  }

  fclose(f);
}

/*
  Function: read_knownGene, similar to read_data but read the format of the
            refFlat file of gene locations
*/
void read_knownGene(const char *filename) {
  FILE * f = fopen(filename, "r");
  if (f == NULL) {printf("Error reading from file %s\n", filename); exit(0); }

  int result=1;
  while (result > 0) {
    int i; char dir;
    RSW_Known * rk = new RSW_Known;
    result = get_line(f, sLine, MAX_LINE);
    if (result < 0) {
      printf("Error reading data file %s, line exceeded %i characters.\n", filename, MAX_LINE);
      delete rk;
      break;
    }
    char *tempA = strtok(sLine, "\t");
    i=0;
    while (tempA != NULL) {
      string temp= tempA; temp.shrink_to_fit();
      pair<unordered_set<string>::iterator,bool> result;
      switch (i) {
      case 0:
        result = stringTable.insert(temp);
        rk->id1 = (result.first)->c_str();
        break;
      case 1:
        result = stringTable.insert(temp);
        rk->id2 = (result.first)->c_str();
        break;
      case 2:
        result = stringTable.insert(temp);
        rk->chromosome = (result.first)->c_str();
        break;
      case 3:
        rk->direction = temp[0];
        break;
      case 4:
        rk->position1 = atol(temp.c_str());
        break;
      case 5:
        rk->position2 = atol(temp.c_str());
        break;
      }
      i++;
      tempA = strtok(NULL,"\t");
      if (tempA == NULL) break;
    }
    if (i < 11) {delete rk; break;}

    if (rk->position1 > rk->position2) {
      int t = rk->position1;
      rk->position1 = rk->position2;
      rk->position2 = t;
    }
    data_known.push_back(*rk);
    delete rk;
  }

  fclose(f);
}

/*
  Function: read_boundaries, similar to read_data but read the format of the
            refFlat file of intron/extron boundaries.
*/
void read_boundaries(const char *filename) {
  FILE * f = fopen(filename, "r");
  if (f == NULL) {printf("Error reading from file %s\n", filename); exit(0); }

  int result=1;
  while (result > 0) {
    int i; char dir;
    RSW_Boundaries * rk = new RSW_Boundaries;
    result = get_line(f, sLine, MAX_LINE);
    if (result < 0) {
      printf("Error reading data file %s, line exceeded %i characters.\n", filename, MAX_LINE);
      delete rk;
      break;
    }
    char *tempA = strtok(sLine, "\t");
    i=0;
    while (tempA != NULL) {
      string temp = tempA; temp.shrink_to_fit();
      pair<unordered_set<string>::iterator,bool> result;
      switch (i) {
        case 0:
          result = stringTable.insert(temp);
          rk->id1 = (result.first)->c_str();
          break;
        case 1:
          result = stringTable.insert(temp);
          rk->id2 = (result.first)->c_str();
          break;
        case 2:
          result = stringTable.insert(temp);
          rk->chromosome = (result.first)->c_str();
          break;
        case 3:
          rk->direction = temp[0];
          break;
        case 11:
          rk->length = atoi(temp.c_str());
          break;
        case 12:
          char *tempB = (char *) malloc(sizeof(char)* (temp.size()+1));
          strcpy(tempB,temp.c_str());
          char *temp1 = strstr(tempB, "--");
          if (temp1 == NULL) {
            rk->position1 = rk->position2 = 0;
          }
          else {
            temp1[0] = '\0';
            rk->position1 = atol(tempB);
            rk->position2 = atol(temp1+2);
          }
          free(tempB);
          break;
      }
      i++;
      tempA = strtok(NULL,"\t");
      if (tempA == NULL) break;
    }
    if (i < 13) {delete rk; break;}

    data_boundaries.push_back(*rk);
    delete rk;
  }

  fclose(f);
}

/*
  Function:   compare_dataById, used for sorting input data

  Sorts based on id, direction, chromosome, position
*/
bool compare_dataById(RSW const &aa, RSW const &bb) {
  int temp = aa.id-bb.id;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  if (aa.direction < bb.direction) 
    return true;
  else if (bb.direction < aa.direction)
    return false;

  temp = aa.chromosome-bb.chromosome;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  if (aa.position < bb.position) return true;
  else return false;
}

/*
  Function:   compare_dataToSort, used for sorting input data

  Sorts based on chromosome, position, id
*/
bool compare_dataByChromPos(RSW const &aa, RSW const &bb) {
  int temp = aa.chromosome - bb.chromosome;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  int position = aa.splitPos - bb.splitPos;
  if (position < 0) return true;
  else if (position > 0) return false;
  
  return false;
}


/*
  Function:  compare_data_known, used for sorting results from refFlat file

  Sort based on chromosome and position.
*/
bool compare_data_known(RSW_Known const &aa, RSW_Known const &bb) {
  int temp = aa.chromosome-bb.chromosome;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  return aa.position1 < bb.position1;
}


/*
  Function: compare_spliceByChromPos, used for sorting junctions

  Sorts based on chromosome, position, splice length - used in sorting
  before computing supporting reads.
*/
bool compare_spliceByChromPos(RSW_splice *aa, RSW_splice *bb) {
  int temp = aa->chromosome-bb->chromosome;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  if (aa->positionSmaller < bb->positionSmaller) return true;
  else if (aa->positionSmaller > bb->positionSmaller) return false;
  
  return false;
}


/*
  Function: readOptionsFromFile, reads parameters for program

  Options assumed to be one per line in a predefined order.
*/
void readOptionsFromFile(const char *filename) {
  FILE * fOptions = fopen(filename, "r");
  if (fOptions == NULL) { printf("Error opening file %s\n", filename); exit(0); }

  int pos = 0; int count = 0;
  char *fields[9]; fields[0] = &options[0];
  int ch;
  while ((ch = fgetc(fOptions)) != EOF) {
    if (pos >= MAX_STR_LEN) { printf("Options file %s is more than the max of %i bytes.\n", filename, MAX_STR_LEN); exit(0); }
    if (ch == '\n') {
      options[pos++] = '\0'; count++;
      if (count < 9)
 fields[count] = &options[pos];
      else break;
    }
    else options[pos++] = ch;
  }

  fclose(fOptions);

  if (count < 9) { printf("Not enough lines in options file.\n"); exit(0); }

  sampleDataFile = fields[0];
  maxDistance = atoi(fields[1]);
  //sampleLength = atoi(fields[2]);
  // this field no longer used, software automatically detects read length for each read in an earlier part of the pipeline
  refFlatFile = fields[3];
  refFlatBoundaryFile = fields[4];
  minSpliceLength = atoi(fields[5]);
  supportPosTolerance = atoi(fields[6]);
  resultsBaseName = fields[7];
  minSupportingReads = atoi(fields[8]);
}

/*
  Function:  setDefaultOptions, sets some default options if now options
             file is given on the command-line.

  Useful for debugging - set the default options to be whatever you are testing,
  so don't have to type in name of options file each time you run the program.
*/
void setDefaultOptions() {
  sampleDataFile = "RSW_test.txt";
  maxDistance = 40000;
  //sampleLength = 33;
  refFlatFile = "refFlat.txt";
  refFlatBoundaryFile = "refFlat.txt.intronBoundary.exonsgaps";
  minSpliceLength = 2;
  supportPosTolerance = 5;
  resultsBaseName = "RSW_tst";
  minSupportingReads = 2;

  printf("Not enough arguments given, using default values.\n");
  printf("Usage is to load options from file: ./splitPairs optionsFile.txt \n");
  printf("And make sure options file has options in order, each on their own line with no extra lines.\n");
}

/*
  Function:  printCurrentOptions, print options to given file pointer.

  Print to file pointer, so can print the options that were used to stdout and/or
  to output files with results.
*/
void printCurrentOptions(FILE *f) {
  fprintf(f, "Running with options...\n"
   "  file with read data          %s\n"
   "  max distance between matches %i\n"
   "  length of samples            variable, auto-detect\n"
   "  refFlat file                 %s\n"
   "  refFlat intron boundary file %s\n"
   "  minimum splice length        %i\n"
   "  tolerance of difference in position for supporting reads  %i\n"
   "  base of file name for writing results                     %s\n"
   "  minimum number of supporting reads                        %i\n\n",
   sampleDataFile, maxDistance, /*sampleLength,*/ refFlatFile, refFlatBoundaryFile, minSpliceLength, supportPosTolerance, resultsBaseName, minSupportingReads);

  if (strlen(sampleDataFile) > MAX_STR_LEN - 100) {
    fprintf(f,"Error, filename %s is too long.\n", sampleDataFile); exit(0);
  }

  fflush(f);
}

// files we will write out to
FILE * fKnown, *fUnknown, //*fKnownFull, *fUnknownFull,
  *fSplitPairs;

// open output files to be ready to write out to them.
void openOutputFiles() {
  sprintf(buff,"%s.results", resultsBaseName);
  fKnown = fopen(buff, "w");
  if (fKnown == NULL) { printf("Error opening file %s for writing.\n", buff); exit(0); }
  printf("Will write summary results that match in known genes to file\n"
  "   %s\n", buff);

  /*sprintf(buff,"%s.results.full", resultsBaseName);
  fKnownFull = fopen(buff, "w");
  if (fKnownFull == NULL) { printf("Error opening file %s for writing.\n", buff); exit(0); }
  printf("Will write full results that match in known genes to file\n"
  "   %s\n", buff);*/

  sprintf(buff,"%s.results.unknown", resultsBaseName);
  fUnknown = fopen(buff, "w");
  if (fUnknown == NULL) { printf("Error opening file %s for writing.\n", buff); exit(0); }
  printf("Will write summary results that do NOT match in known genes to file\n"
  "   %s\n", buff);

  /*  sprintf(buff,"%s.results.unknown.full", resultsBaseName);
  fUnknownFull = fopen(buff, "w");
  if (fUnknownFull == NULL) { printf("Error opening file %s for writing.\n", buff); exit(0); }
  printf("Will write full results that do NOT match in known genes to file\n"
  "   %s\n", buff);*/

  sprintf(buff,"%s.results.splitPairs", resultsBaseName);
  fSplitPairs = fopen(buff, "w");
  if (fSplitPairs == NULL) { printf("Error opening file %s for writing.\n", buff); exit(0); }
  printf("Will write split pairs to file\n"
  "   %s\n", buff);
}

/*
  Function:  printSplice, print a given junction to the given opened file
 Commenting this out and printing out one that is more comparable to the old output
   - Aaron
*/
/*void printSplice(FILE *f, RSW_splice *sp) {
  fprintf(f,
      "%s\t%s\t%li\t%li\t%li\t%li\t%li--%li\t%li--%li\t%s", 
      sp->geneName, sp->chromosome,
      sp->numSupport,
      sp->numSupportHalves,
      sp->numSupportTotal,
      sp->positionLarger-sp->positionSmaller,
      sp->positionSmaller,sp->positionLarger,
      sp->minSmallSupport,sp->maxLargeSupport,
      sp->novel ? "Novel" : "*"
  );
}
*/

void printSplice(FILE *f, RSW_splice *sp) {
  fprintf(f,
      "%s\t%s\t%li\t%li\t%li\t%li\t%li--%li\t%s", 
      sp->geneName, 
      sp->chromosome,
      sp->numSupport, //full support 
      sp->numSupportHalves, //half support
      sp->numSupportTotal, //"total" support
      sp->positionLarger - sp->positionSmaller, //splice length
      sp->minSmallSupport,sp->maxLargeSupport, //range of indices jnctn occrs
      sp->novel ? "Novel" : "*"
  );
}



string getHalfStats() {
  char s[10000];
  
  int halfCount = 0;
  int minMax = -1, minMin = -1, maxMax = -1, maxMin = -1, minTotal = 0, maxTotal = 0;
  for(auto it=data_halves.begin(); it != data_halves.end(); it++) {
    halfCount++;
    int min = it->second.minLength, max = it->second.maxLength;
    minTotal += min;  maxTotal += max;
    if (minMax == -1 || min > minMax) minMax = min;
    if (minMin == -1 || min < minMin) minMin = min;
    if (maxMax == -1 || max > maxMax) maxMax = max;
    if (maxMin == -1 || max < maxMin) maxMin = max;
  }
  sprintf(s, "Half lengths:                               min %lf avg, range %li-%li;  max %lf avg, range %li-%li",
   (double) minTotal / halfCount, minMin, minMax,
   (double) maxTotal / halfCount, maxMin, maxMax);

  return s;
}



/*
  Function:  printStats, prints statistics gathered so far to the
             opened file - useful for debugging to see some partial
             information as each phase of the program finishes.
*/
void printStats(FILE *f) {
  FILE *fStatus = fopen("/proc/self/status","r");
  char s[1000], mem[20]="", units[20]="";
  while (fscanf(fStatus,"%999s",s) == 1) {
    if (strcmp(s,"VmRSS:") == 0) {
      int result = fscanf(fStatus,"%19s %19s",mem, units);
      break;
    }
  }
  fclose(fStatus);

  endTime = time(NULL);
  fprintf(f, "Finished processing data, results written to files.\n");
  fprintf(f, "Number of entries in data file:             %li\n", data.size());
  fprintf(f, "Number of different reads:                  %i\n", numDifferentReads);
  fprintf(f, "Number of entries in refFlat file:          %li\n", data_known.size());
  fprintf(f, "Number of entries in refFlat boundary file: %li\n", data_boundaries.size());
  fprintf(f, "Number of matches:                          %li\n", data_splice.size());
  fprintf(f, "String table size:                          %li\n", stringTable.size());
  fprintf(f, "VmRSS, memory resident set size:            %s %s\n", mem, units);
  fprintf(f, "Total time to process:                      %li seconds\n", endTime-beginTime);

  // statistics for halves...
  fprintf(f, "%s", halfStatsString.c_str());
  fprintf(f, "\n");

  fflush(f); // force write to disk
}


/*
  Check a half in data against a splice in data_splice, in particular
  check data_splice[sp1] against data[i_data].
  Return: 1 if they are a match
          0 if not a match
   -1 if should break out of loop back in main (stop incrementing i_data because past data_splice[sp1] in data)
 */
int checkHalf(int sp1, int i_data, bool smallEnd) {
  int c = data_splice[sp1]->chromosome - data[i_data].chromosome;
  
  // if not same chromosome, either wait for sp1 to catch up, or let i_data catch up
  if (c < 0) return -1;
  else if (c == 0) {
    int p;
    if (smallEnd) p = data_splice[sp1]->positionSmaller - data[i_data].splitPos;
    else p = data_splice[sp1]->positionLarger - data[i_data].splitPos;
    if ( p < 0) return -1;
    else if (p == 0) { // a match
      auto fOther = data_halves.find(data[i_data].otherHalfKey);
      if (data_splice[sp1]->direction == data[i_data].direction && 
   (fOther == data_halves.end() ||
    fOther->second.maxLength < data[i_data].totalReadLength - data[i_data].length)) {
 // then this is a half that is at the right position and doesn't have a matching
 // other half (presumably because of being in the max file) that is long enough so let's count it.
 return 1;
      }
    }
    else // p > 0
      return 0;
  }
  else // c > 0
    return 0;
}


int main(int argc, char *argv[]) {
  setpriority(0, 0, 20); // so other processes get priority over this one

  beginTime = time(NULL);

  // read options, from file or default options
  if (argc > 1) 
    readOptionsFromFile(argv[1]);
  else 
    setDefaultOptions();

  // write out options to all output files and stdout
  openOutputFiles();

  printCurrentOptions(stdout);
  printCurrentOptions(fKnown);
  //  printCurrentOptions(fKnownFull);
  printCurrentOptions(fUnknown);
  //  printCurrentOptions(fUnknownFull);
  printCurrentOptions(fSplitPairs);


  // read from refFlat file into data_known array, 
  read_knownGene(refFlatFile);
  sort(data_known.begin(), data_known.end(), compare_data_known);
  printf("Done reading/sorting refFlat, total time elapsed %li seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // read from refFlat boundary file into data_boundaries array, 
  read_boundaries(refFlatBoundaryFile);
  printf("Done reading refFlat intron/exon boundaries, total time elapsed %li seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // read the read data
  read_data(sampleDataFile);
  halfStatsString = getHalfStats();
  printf("Done reading read data, total time elapsed %li seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // sort the read data
  sort(data.begin(), data.end(), compare_dataById);
  printf("Done sorting read data, total time elapsed %li seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // look at all pairs of read segments, looking for matches
  numDifferentReads = 0;
  const int inputSize = data.size();
  for(int left=0; left < inputSize; left++) {
    int right;
    if (left == 0 || (data[left].id  != data[left-1].id) ) {
      numDifferentReads++;
    }

    for(right=left+1; right < inputSize; right++) {
      // read segments are ordered by id/chromosome/strand, so if 
      // there isn't a match we can skip the rest of the read segments
      // for "right", and go to the next iteration of the "left" loop
      if (data[left].id != data[right].id ||
          data[left].direction != data[right].direction ||
          data[left].chromosome != data[right].chromosome) {
        break;
      }

      if (data[right].position - data[left].position > maxDistance) {
        break;
      }

      // want it to be from two sides of the same segment
      if (data[left].side == data[right].side) {
        continue;
      }

      // total length should be correct
      if (data[left].length + data[right].length != data[left].totalReadLength) {
        continue;
      }

      // calculate the end of the segments, since what is given
      // in the data is the beginning of the segments
      int endSmaller, endLarger; // splice is between endSmaller and endLarger
      int first, second;
      if (data[left].side == 'L' && data[left].direction == '+' ||
          data[left].side == 'R' && data[left].direction == '-') { 
        first = left; second = right;
      }
      else { 
        first = right; second = left;
      }
      endSmaller = data[first].position + data[first].length;
      endLarger = data[second].position;

      // splice length, and check that it is within specified bounds
      int spliceLength = endLarger - endSmaller;
      if (spliceLength > maxDistance) continue;
      if (spliceLength < minSpliceLength) continue;

      // check if we already have this splice from this read...
      int i;
      for(i=0; i < data_splice.size(); i++) {
        if (data_splice[i]->positionSmaller != endSmaller ||
          data_splice[i]->positionLarger != endLarger ||
          data_splice[i]->id != data[left].id ||
          data_splice[i]->chromosome != data[left].chromosome)
        {
          continue;
        }
        break;
      }
      // if already have this exact splice for this chromosome from this read, don't include it again.
      if (i < data_splice.size()) continue;

      // note: could print this match here, step 5 done.

      // look for this in the known gene...
      int k; int foundInGene = 0;
      for(k=0; k < data_known.size(); k++) {
        if ((((data[left].chromosome == data_known[k].chromosome) &&
          (data[left].position >= data_known[k].position1 &&
           data[right].position >= data_known[k].position1) &&
          (data[left].position <= data_known[k].position2 &&
           data[right].position <= data_known[k].position2)))) 
        {
          if (foundInGene) ;//printf("DUPLICATE_"); // duplicate //TODO: remove this? -aaron
          foundInGene = 1;
          break; // just cut off search, don't look for duplicates
        }
      }
      int geneIndex = k;
      
      // make a new splice record and put into vector of splices
      RSW_splice *sp = new RSW_splice;
      if (sp == NULL) {
        printf("ERROR, new in C++ failed, maybe out of memory.\n");
        exit(0);
      }
      sp;
      if (foundInGene) {
        sp->geneName = data_known[geneIndex].id1;
        sp->geneUnknown = 0;
      } 
      else {
        sp->geneName = unfound_string;
        sp->geneUnknown = 1;
      }
      sp->id = data[left].id;
      sp->chromosome = data[left].chromosome;
      sp->direction = data[left].direction;
      sp->positionSmaller = sp->minSmallSupport = endSmaller;
      sp->positionLarger = sp->maxLargeSupport = endLarger;
      sp->alreadyReported = false;
      sp->print = false;
      sp->numSupport = sp->numSupportHalves = sp->numSupportTotal = 0;
      sp->leftLength = data[left].length;

      data_splice.push_back(sp);
    }
  }

  printf("Done finding matched pairs, total time elapsed %li seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // re-sort input data by chromosome and position
  sort(data.begin(), data.end(), compare_dataByChromPos);
  printf("Done resorting input data by chromosome and position, total time elapsed %li seconds\n", time(NULL)-beginTime);
  printStats(stdout);
  
  
  // sort splices by chromosome and position
  sort(data_splice.begin(), data_splice.end(), compare_spliceByChromPos);

  printf("Done sorting matched pairs, total time elapsed %li seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  int i_lastEndSmaller = 0, i_lastEndLarger = 0;

  // compute supporting reads, print them out on the fly.

  for(int sp1=0; sp1 < data_splice.size(); sp1++) {
    int sp2;

    if (data_splice[sp1]->alreadyReported) continue;
    
    unordered_set <const char *> supported_read_ids; // list of supporting reads
    unordered_set <const char *> supported_read_ids_halves; // list of supporting reads
    unordered_set <const char *> supported_read_ids_both;   // list of supporting reads
    unordered_set <RSW_splice *> supported_splices;
    //unordered_set <int> supported_halves; // int is the index into data - note that only works as long as data is not resorted
    supported_read_ids.insert(data_splice[sp1]->id);
    supported_read_ids_both.insert(data_splice[sp1]->id);
    supported_splices.insert(data_splice[sp1]);

    // scan through the following reads ...
    for(sp2=sp1; sp2 < data_splice.size(); sp2++) {
      // see if these reads support each other

      // if not same chromosome, no support.
      if (data_splice[sp1]->chromosome != data_splice[sp2]->chromosome) {
        break;
      }

      // only need to go up to supportPosTolerance away in position, then break
      if (data_splice[sp2]->positionSmaller > data_splice[sp1]->positionSmaller + supportPosTolerance) {
        break;
      }

      if (abs(data_splice[sp1]->positionLarger - data_splice[sp1]->positionSmaller) != 
          abs(data_splice[sp2]->positionLarger - data_splice[sp2]->positionSmaller)) {
        continue; // splice length must be the same
      }
      
      // note: if id already 
      if (data_splice[sp2]->id == data_splice[sp1]->id) {
        // if report sp1, then shouldn't report other splices for sp1 that are close
        supported_splices.insert(data_splice[sp2]);
        continue; 
      }

      // they are matches for each other
      supported_read_ids.insert(data_splice[sp2]->id); 
      supported_read_ids_both.insert(data_splice[sp2]->id); 
      supported_splices.insert(data_splice[sp2]);
      
      if (data_splice[sp2]->positionSmaller < data_splice[sp1]->minSmallSupport) {
        data_splice[sp1]->minSmallSupport = data_splice[sp2]->positionSmaller;
      }

      if (data_splice[sp2]->positionLarger > data_splice[sp1]->maxLargeSupport) {
        data_splice[sp1]->maxLargeSupport = data_splice[sp2]->positionLarger;
      }
    }

    // scan also through read data, looking for halves that
    // match up but don't have the other half because it is probably
    // in the max file.
    for(; ; i_lastEndSmaller++)  {
      int result = checkHalf(sp1, i_lastEndSmaller, true);
      if (result < 0) {
        break;
      }
      else if (result > 0) {
        // then this is a half that is at the right position and doesn't have a matching
        // other half (presumably because of being int he max file), so let's count it.
        supported_read_ids_halves.insert(data[i_lastEndSmaller].id); 
        supported_read_ids_both.insert(data[i_lastEndSmaller].id);
      }
    }

    // and similarly, check the larger half of the split.
    for(int i_data = i_lastEndSmaller; ; i_data++)  {
      int result = checkHalf(sp1, i_data, false);
      if (result < 0) break;
      else if (result > 0) {
        // then this is a half that is at the right position and doesn't have a matching
        // other half (presumably because of being int he max file), so let's count it.
        supported_read_ids_halves.insert(data[i_data].id); 
        supported_read_ids_both.insert(data[i_data].id);
      }
    }

    
    if (supported_read_ids_both.size() >= minSupportingReads) {
      data_splice[sp1]->print = true;
      data_splice[sp1]->alreadyReported = true;
      data_splice[sp1]->numSupport = supported_read_ids.size();
      data_splice[sp1]->numSupportHalves = supported_read_ids_halves.size();
      data_splice[sp1]->numSupportTotal = supported_read_ids_both.size();

      for(const auto& x: supported_splices) {
 x->alreadyReported = true;
      }       
    }
  }

  printf("Done computing supporting reads, total time elapsed %li seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // sort splices by chromosome and length, and # supporting reads
  //sort(data_splice.begin(), data_splice.end(), compare_spliceByChromLen);
  //printf("Done sorting matched pairs again, total time elapsed %li seconds\n", time(NULL)-beginTime);
  //printStats(stdout);

  //printf("Done filtering matched pairs, total time elapsed %li seconds\n", time(NULL)-beginTime);
  //printStats(stdout);

  endTime = time(NULL);


  printf("Done calculating/filtering supporting reads, total time elapsed %li seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // print the splice results.  start with the statistics
  printStats(stdout);
  printStats(fKnown);
  //  printStats(fKnownFull);
  printStats(fUnknown);
  //  printStats(fUnknownFull);
  printStats(fSplitPairs);

  int k;


  // save all the splices, and reads that support each one into .splitPairs file - this is
  // the full results, with many duplicates of splices.  The file also can be HUGE, so
  // generally this file will often get deleted unless needed for debugging.
  fprintf(fSplitPairs, "Id\tGene\tChr\t# Supporting reads\t# Supporting halves\t# Supporting total\tLength\tSplice region\tSupporting splice range\tLeft side length\n");
  int i_data=0;
  for(k=0; k < data_splice.size() && i_data < data.size();) {
    if (i_data == data.size() ||
 data_splice[k]->chromosome < data[i_data].chromosome ||
 (data_splice[k]->chromosome == data[i_data].chromosome && 
  data_splice[k]->positionSmaller < data[i_data].splitPos)) {
      // print a splice
      fprintf(fSplitPairs, "%s\t%s\t%s\t%li\t%li\t%li\t%li\t%li-%li\t%li-%li\t%li\n", 
       data_splice[k]->id, data_splice[k]->geneName,
       data_splice[k]->chromosome,
       data_splice[k]->numSupport,
       data_splice[k]->numSupportHalves,
       data_splice[k]->numSupportTotal,
       data_splice[k]->positionLarger-data_splice[k]->positionSmaller,
       data_splice[k]->positionSmaller,data_splice[k]->positionLarger,
       data_splice[k]->minSmallSupport,data_splice[k]->maxLargeSupport,
       data_splice[k]->leftLength
       );
      k++;
    }
    else {
      // print a half that doesn't have a matching other half that is big enough
      auto fOther = data_halves.find(data[i_data].otherHalfKey);
      if (fOther == data_halves.end() ||
          fOther->second.maxLength < data[i_data].totalReadLength - data[i_data].length) {
        fprintf(fSplitPairs, "%s\t%s\t%s\t%li\t%li\t%li\t%li\t%li-%li\t%li-%li\t%li %c %c\n",
          data[i_data].id, "???",
          data[i_data].chromosome,
          0,0,0,
          0,
          data[i_data].position,data[i_data].splitPos,
          0,0,
          data[i_data].length,
          data[i_data].side, data[i_data].direction
        );
      }
      i_data++;
    }
  }

  // save the tabulated results.  the fKnown file is the only one normally looked at.
  // fKnownFull has the same information as fKnown, but also has the list of supporting reads for each junction.
  // fUnknown and fUnknownFull are for junctions not within genes.
  fprintf(fKnown, "GeneName\tChromosome\t# supporting reads\t# supporting halves\t# supporting total\tsplice length\trange of supporting reads\tNovel or not (*)\n"); 
  fprintf(fUnknown, "GeneName\tChromosome\t# supporting reads\t# supporting halves\t# supporting total\tsplice length\trange of supporting reads\tNovel or not (*)\n"); 
  FILE * f;//, *fFull;
  for(k=0; k < data_splice.size(); k++) {
    // already decided if we should print this one or not
    if (! data_splice[k]->print) continue;

    // check if novel or not.
    int j;
    data_splice[k]->novel = true;
    for(j=0; j < data_boundaries.size(); j++) {
      if (data_splice[k]->positionLarger-data_splice[k]->positionSmaller != 
          data_boundaries[j].length) {
        continue;
      }
      if (abs(data_splice[k]->minSmallSupport-data_boundaries[j].position1) <= supportPosTolerance &&
          abs(data_splice[k]->maxLargeSupport-data_boundaries[j].position2) <= supportPosTolerance) {
        break;
      }
    }
    if (j < data_boundaries.size()) { // if found in the boundaries data, not novel
      data_splice[k]->novel = false;
    }

    // going into known file or unknown
    if (data_splice[k]->geneUnknown) {
      f = fUnknown; //fFull = fUnknownFull;
    }
    else {
      f = fKnown; //fFull = fKnownFull;
    }

    // print out results
    printSplice(f, data_splice[k]);
    fprintf(f, "\n");
  }

  fclose(fKnown); //fclose(fKnownFull);
  fclose(fUnknown); //fclose(fUnknownFull);
  fclose(fSplitPairs);

  printf("Done saving results, total time elapsed %li seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // free memory.  good to do so we can run a memory checker and verify
  // we don't have any memory leaks.
  while (data_splice.size() > 0) {
    delete data_splice.back(); data_splice.pop_back();
  }
  
  return 0;
}
