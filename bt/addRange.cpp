
#define _XOPEN_SOURCE 700

#include "addRange.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <string>
#include <unistd.h>
#include <ctype.h>
#include <strings.h>


const std::string TMP_FILENAME = ".w_SpliceSequence";
const std::string HEADER_ID = "genename";
/*const std::string RANGE_COL = "range of supporting reads";
const std::string BRACKETSEQ_COL = "Bracketed sequence";
const std::string CHR_COL = "Chromosome";
const std::string SPLICELEN_COL = "splice length";*/
const int BOUNDARY_LEN = 15;


const int GENE_COL = 0;
const int CHR_COL = 1;
const int SUPPORT_READS_COL = 2;
const int SPLICE_LEN_COL = 3;
const int RANGE_COL = 4;
const int NOVEL_COL = 5;
const int BRACKET_SEQ_COL = 6;
const int SPLICED_SEQ_COL = 7;
 

void usage();
void bracketRange(std::string &range);
void modifyResults(std::ifstream& resultFile, std::ofstream& writeFile, std::string genFileName);
void extendRange(std::string token, std::string &bracketL, std::string &bracketR, long &moveLeft, long &moveRight);
//int getColumnNo( std::string header, std::string column );
void makeNewStrs( std::string &newBrackSeq, std::string &newSpliceSeq, FILE *spliceFile, int spliceLen, int boundaryLen );
/*void makeNewStrs_saveEdges( std::string &newBrackSeq, std::string &newSpliceSeq, FILE *spliceFile, int bp, long wideL,
                  long rangeL, long rangeR, long wideR, std::string &bracketL, std::string &bracketR );*/
std::string seekHeaderLine( std::ifstream& resultFile, std::ofstream& writeFile );
std::string editRange(std::string& range);
std::vector<std::string> splitString(std::string str, const char delimit);



/* Prints brief usage message and exits the program. */
void usage() {
	std::cout << "usage: addRange resultsFile genomeFile";
	std::cout << std::endl << std::endl;
	exit(EXIT_FAILURE);
	return;
}

/* Rough equivalent of strtok.  Strips away newline ('\n') characters as well.
 * Returns the tokens in a vector<string> */
std::vector<std::string> splitString( std::string str, const char delimit ) {
	std::vector<std::string> tokens;
	std::string buf;
	for( int i = 0; i < str.length(); i++ ) {
		if( str[i] != delimit && str[i] != '\n') {
			buf.append(1, str[i]);	
		}
		else {
			tokens.push_back(buf);
			buf = std::string();
		}
	}
	if( !buf.empty() ) tokens.push_back(buf);
	return tokens;
}

/* Finds the column number associated with a substring of a column header
 * from the argument header string.  Returns -1 if range column was not found.
 * Otherwise, returns the range of supporting reads column. 
 * WARNING: This function is currently case sensitive when comparing
 * TO DO: make this fx NOT case sensitive. */
/*int getColumnNo( std::string header, std::string column ) {
	int currCol = 0;
	int rangeCol = -1;

	std::vector<std::string> tokens = splitString(header, '\t');
	for( std::vector<std::string>::iterator itr = tokens.begin(); itr != tokens.end(); ++itr ) {
		if( itr->find(column) != std::string::npos ) {
			if( rangeCol == -1 ) {
				rangeCol = currCol;
			}
		}
		currCol++;
	}
	return rangeCol;
}*/

/* This is a little flimsy: Find the header based on the first column: GeneName
 * Return the header line string, with all extra columns appended. 
 * Return empty string if header line wasn't found. */
std::string seekHeaderLine( std::ifstream &resultFile, std::ofstream &writeFile ) {
	bool headerFound = false;
	std::string buf;

	while( !headerFound ) {
		// Grab the next line from the input result file
		std::getline(resultFile, buf);
		if( resultFile.eof() ) {
			break;
		}
		else {
			// Check to see if the line is the header line.
			if( strncasecmp(buf.c_str(),HEADER_ID.c_str(),HEADER_ID.length()) == 0 ) {
				// If it is, append all extra header values
				buf = buf.append("\tsplice sequence");
				headerFound = true;
			}
			// Print the file to writefile
			writeFile << buf << std::endl;
			writeFile.flush();
		}
	}
	
	if( headerFound ) {
		return buf;
	}
	return std::string();  //else return empty string
}

/* Takes the string &range and edits it from ABCD--EFGH to ABCD]--[EFGH form.
 * If the form isn't in the correct form to begin with, this function does
 * nothing and returns immediately. */
void bracketRange(std::string &range) {
	std::string newRange;
	int left = range.find("--");
	if( left == std::string::npos ) return;
	int right = left + 2;
	newRange = newRange.append(range.substr(0,left)).append("]--[").append(range.substr(right));
	range = newRange;
	return;
}
/* Edit the first number in "range of supporting reads" and add brackets
 * to show bracketed inclusion decisively for easier reading */
std::string editRange(std::string &range) {
	bracketRange(range);
	std::string newRange;
	int endFirstNo = range.find("]--[");
	newRange = newRange.append( std::to_string( atol(range.substr(0,endFirstNo).c_str()) - 1 ) );
	newRange = newRange.append( range.substr(endFirstNo) );
	return newRange;
}


/* All arguments are references but token.  This function takes the bracketed sequence
 * string and splits the left and right sides into bracketL and bracketR.
 * The sizes of these left and right strings are stored in moveLeft and moveRight */
void extendRange(std::string token, std::string &bracketL, std::string &bracketR, long &moveLeft, long &moveRight) {
	int pos = token.find("--");
	bracketL = token.substr(0,pos);
	bracketR = token.substr(pos+2);
	moveLeft = bracketL.length();
	moveRight = bracketR.length();
	return;
}

/* ALERT:  the current preferred method is to take an arbitrary number of outside
 * nucleotides, so this (mostly accurate) method shouldn't be used unless that changes
 *
 * The FILE* spliceFile gets its input from bowtie-inspect-RSR command, which takes
 * a maximal possible string comprising of the range of supporting reads with the 
 * size of each side of the bracketed sequence added to either side of the range.
 * 
 * As input, this function takes the original range values, the maximally wide values,
 * the length of the supported range, and removes the overlap the given bracketed sequence
 * has with the range of supporting range.  If the (splice length == range.length), then
 * there will be nothing removed from the bracketed sequence.
 *
 * There is no return value, but the new bracketed sequence is stored in newBrackSeq and
 * the full spliced sequence (range) is stored in newSpliceSeq.  */
/*void makeNewStrs_saveEdges( std::string &newBrackSeq, std::string &newSpliceSeq, FILE *spliceFile, int spliceLen, long wideL,
                  long rangeL, long rangeR, long wideR, std::string &bracketL, std::string &bracketR ) {
	// First, read all input from spliceFile and store it into a spliceStr string
	std::string spliceStr;
	char spliceBuf[wideR - wideL + 2];
	char c;
	int i = 0;
	while( (c = fgetc(spliceFile)) != EOF ) {
		if( c == 'A' || c == 'C' || c == 'G' || c == 'T' ) {
			spliceBuf[i++] = c;
		}
	}
	spliceBuf[i] = '\0';
	spliceStr = std::string(spliceBuf);
	// Now, pull out the left- and right-bracketed read regions
	std::size_t lpos;	
	std::string leftStr, rightStr;
	lpos = spliceStr.find(bracketL);
	leftStr = spliceStr.substr(lpos, bracketL.length() - lpos);
	spliceStr = spliceStr.substr(bracketL.length());
	rightStr = spliceStr.substr(rangeR - rangeL);
	spliceStr = spliceStr.substr(0, rangeR - rangeL);
	rightStr = rightStr.substr(0, rightStr.length() - (rangeR - rangeL - spliceLen - lpos));

	newBrackSeq = leftStr.append("]--[").append(rightStr);
	newSpliceSeq = spliceStr;
	return;
}*/

void makeNewStrs( std::string &newBrackSeq, std::string &newSpliceSeq, FILE *spliceFile, int spliceRangeLen, int boundaryLen ) {
	// First, read all input from spliceFile and store it into a spliceStr string
	std::string spliceStr, leftStr, rightStr;
	char spliceBuf[spliceRangeLen + 2*boundaryLen + 2];
	char c;
	int i = 0;
	while( (c = fgetc(spliceFile)) != EOF ) {
		if( c == 'A' || c == 'C' || c == 'G' || c == 'T' ) {
			spliceBuf[i++] = c;
		}
	}
	spliceBuf[i] = '\0';
	spliceStr = std::string(spliceBuf);
	leftStr = spliceStr.substr(0,boundaryLen);
	rightStr = spliceStr.substr( spliceStr.length() - boundaryLen );
	spliceStr = spliceStr.substr( boundaryLen, spliceRangeLen );

	newBrackSeq = leftStr.append("]--[").append(rightStr);
	newSpliceSeq = spliceStr;
}

/* This is the part of the file that is changing the results file from the original
 * sp4 form to the version with the Spliced Sequence. */
void modifyResults(std::ifstream &resultFile, std::ofstream &writeFile, std::string genFileName) {
	std::string headerStr = seekHeaderLine(resultFile, writeFile);	
	std::vector<std::string> tokens;
	std::string buf;

	getline(resultFile,buf);
	while( !resultFile.eof() ) {
		int spliceLen;
		long wideL, wideR, rangeL, rangeR;
		tokens = splitString( buf, '\t' );
		std::string chr, bracketL, bracketR;
		
		int tokNo = 0;
		// In this iteration through a line's data sections, grab information but do not 
		// print or edit contents of the tokens
		for(std::vector<std::string>::iterator itr = tokens.begin(); itr != tokens.end(); ++itr ) {
			if( tokNo == CHR_COL ) {
				chr = std::string(*itr);
			}
			if( tokNo == BRACKET_SEQ_COL ) {
				extendRange(*itr, bracketL, bracketR, wideL, wideR);
			}
			if( tokNo == RANGE_COL ) {
				int pos = itr->find("--");
				rangeL = atol( itr->substr(0,pos).c_str() );
				rangeR = atol( itr->substr(pos+2).c_str() );
			}
			if( tokNo == SPLICE_LEN_COL ) {
				spliceLen = atol( itr->c_str() );
			}
			tokNo += 1;
		}
		// Set wideL & wideR so they define position explicitly instead of a move-by value
		wideL = rangeL - BOUNDARY_LEN;
		wideR = rangeR + BOUNDARY_LEN;

		// build inspect-RSR command
		std::string rsrCmd = std::string("bowtie-inspect-RSR -c ");
		rsrCmd = rsrCmd.append(chr).append(" -i ").append( std::to_string(wideL) );
		rsrCmd = rsrCmd.append(" -j ").append( std::to_string(wideR) ).append(" ");
		rsrCmd = rsrCmd.append(genFileName);
		
		//grab string from bowtie-inspect-RSR
		FILE* spliceFile = popen(rsrCmd.c_str(), "r");
		if( spliceFile == NULL ) {
			std::cerr << "Unable to run bowtie-inspect-RSR" << std::endl;
			exit(EXIT_FAILURE);
		}

		std::string newBrackSeq, newSpliceSeq;
		makeNewStrs(newBrackSeq, newSpliceSeq, spliceFile, rangeR - rangeL, BOUNDARY_LEN );
		tokNo = 0;
		//print the modified current line to the writefile
		for(std::vector<std::string>::iterator itr = tokens.begin(); itr != tokens.end(); ++itr ) {
			if( tokNo == RANGE_COL ) {
				// Edit the range of supp. reads column to show proper inclusion w/ square brackets
				*itr = editRange(*itr);
			}
			if( tokNo == BRACKET_SEQ_COL ) {
				*itr = newBrackSeq;
			}
			writeFile << *itr << '\t';
			writeFile.flush();
			tokNo += 1;
		}
		// Write the spliced sequence to the file and grab the next line
		writeFile << newSpliceSeq << std::endl;
		writeFile.flush();
		pclose(spliceFile);
		getline(resultFile,buf);
	} 
	return;
}

int main(int argc, char *argv[]) {
	int currArg;
	int rangeCol = -1;
	std::string resultFileName;
	std::string genFileName;
	std::string writeFileName;
	std::ifstream resultFile;
	std::ofstream writeFile;
	
	while( (currArg = getopt(argc, argv, "c:")) != -1 ) {
		switch(currArg) {
			case 'c':
				/*rangeCol = atoi(optarg);*/
				break;
			default:
				usage();
			}
	}
	// Process results and genome files.  There should only be two files remaining.
	if( argc - 2 != optind ) usage();
	currArg = optind;
	resultFileName = std::string(argv[currArg]);
	resultFile.open(resultFileName, std::ifstream::in);
	if( !resultFile.good() ) {
		std::cerr << "Error opening result file." << std::endl;
		exit(EXIT_FAILURE);
	}
	genFileName = std::string(argv[++currArg]);
	writeFileName = std::string(resultFileName).append(TMP_FILENAME);
	writeFile.open(writeFileName, std::ofstream::out);
	if( !writeFile.good() ) {
		resultFile.close();
		std::cerr << "Error opening output file." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	modifyResults(resultFile, writeFile, genFileName);
	resultFile.close();
	writeFile.close();
	return 0;
}
