/* This file is used to create bowtie-inspect-RSR, which extends the functionality of
 * the bowtie-inspect program to allow for the extraction of specific assembly sequences
 * and for the processing of RSR/RSF-specific output files.
 *
 * Contributing authors to additions and modifications in this file:
 *
 * Jeff Kinne <jkinne@cs.indstate.edu>
 * Aaron Cox  <acox@cs.indstate.edu>
 */


#include <string>
#include <vector>
#include <iostream>
#include <getopt.h>
#include <fstream>
#include <stdexcept>
#include <seqan/find.h>

#include "assert_helpers.h"
#include "endian_swap.h"
#include "ebwt.h"
#include "reference.h"
#include "addRange.h"

using namespace std;
using namespace seqan;

static bool showVersion = false; // just print version and quit?
static bool batchMode = false; // If true, run ISU-specific version of -c/-i/-j mode on a file
int verbose             = 0;  // be talkative
static int names_only   = 0;  // just print the sequence names in the index
static int summarize_only = 0; // just print summary of index and quit
static int across       = 60; // number of characters across in FASTA output
static bool extra       = false; // print extra summary info
static bool exclAllGaps = false; // print extra summary info
static bool refFromEbwt = false; // true -> when printing reference, decode it from Ebwt instead of reading it from BitPairReference
static long startIndex = 0;
static long stopIndex = -1;
static string chrName = "";
static string resultFileName = "";
static string outputFileName = "";
static string wrapper;
static const char *short_options = "vhnsea:i:j:c:f:o:";

enum {
	ARG_VERSION = 256,
	ARG_USAGE,
	ARG_EXTRA,
	ARG_EXCL_AMBIG,
	ARG_WRAPPER
};

static struct option long_options[] = {
	{(char*)"verbose",  no_argument,        0, 'v'},
	{(char*)"version",  no_argument,        0, ARG_VERSION},
	{(char*)"usage",    no_argument,        0, ARG_USAGE},
	{(char*)"extra",    no_argument,        0, ARG_EXTRA},
	{(char*)"excl-ambig",no_argument,       0, ARG_EXCL_AMBIG},
	{(char*)"names",    no_argument,        0, 'n'},
	{(char*)"summary",  no_argument,        0, 's'},
	{(char*)"help",     no_argument,        0, 'h'},
	{(char*)"across",   required_argument,  0, 'a'},
	{(char*)"ebwt-ref", no_argument,        0, 'e'},
	
	{(char*)"start-index", required_argument,     0, 'i'},
	{(char*)"stop-index", required_argument,      0, 'j'},
	{(char*)"chr-name", required_argument,        0, 'c'},
	{(char*)"result-file", required_argument,        0, 'f'},
	{(char*)"output-file", required_argument,        0, 'o'},
	
	{(char*)"wrapper",  required_argument,  0, ARG_WRAPPER},
	{(char*)0, 0, 0, 0} // terminator
};

/**
 * Print a summary usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out
	  << "Note - modified bowtie at Indiana State University..." << endl
	<< "Usage: bowtie-inspect [options]* <ebwt_base>" << endl
	<< "  <ebwt_base>        ebwt filename minus trailing .1." + gEbwt_ext + "/.2." + gEbwt_ext << endl
	<< endl
	<< "  By default, prints FASTA records of the indexed nucleotide sequences to" << endl
	<< "  standard out.  With -n, just prints names.  With -s, just prints a summary of" << endl
	<< "  the index parameters and sequences.  With -e, preserves colors if applicable." << endl
	<< endl
	<< "Options:" << endl;
	if(wrapper == "basic-0") {
		out << "  --large-index      force inspection of the 'large' index, even if a" << endl
			<< "                     'small' one is present." << endl;
	}
	out
	  << "  -f/--result-file <text>  File containing entries with chromosome, start-index, stop-index" << endl
	  << "                             Running -f option with no -o option outputs to stdout." << endl
	  << "  -o/--output-file <text>  File to output modified version of result-file with spliced sequence to." << endl
	  << "                             Specify 'default' to append '.splSeq' to end of result-file name." << endl
	  << "  -i/--start-index <int>   Index in chromosome(s) to start printing from (default: 0)" << endl
	  << "  -j/--stop-index <int>    Index in chromosome(s) to stop printing at (default: -1)" << endl
	  << "  -c/--chr-name <text>     Name of chromosome to print from (default: blank)" << endl
	  << "  -a/--across <int>        Number of characters across in FASTA output (default: 60)" << endl
		<< "  -n/--names               Print reference sequence names only" << endl
		<< "  -s/--summary             Print summary incl. ref names, lengths, index properties" << endl
		<< "  -e/--ebwt-ref            Reconstruct reference from ebwt (slow, preserves colors)" << endl
		<< "  -v/--verbose             Verbose output (for debugging)" << endl
		<< "  -h/--help                print detailed description of tool and its options" << endl
		<< "  --help                   print this usage message" << endl
	;
	if(wrapper.empty()) {
		cerr << endl
		     << "*** Warning ***" << endl
			 << "'boowtie-inspect' was run directly.  It is recommended "
			 << "to use the wrapper script instead."
			 << endl << endl;
	}
}

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
static int parseInt(int lower, const char *errmsg) {
	long l;
	char *endPtr= NULL;
	l = strtol(optarg, &endPtr, 10);
	if (endPtr != NULL) {
		if (l < lower) {
			cerr << errmsg << endl;
			printUsage(cerr);
			throw 1;
		}
		return (int32_t)l;
	}
	cerr << errmsg << endl;
	printUsage(cerr);
	throw 1;
	return -1;
}

string openFile(char *fileName, fstream file, char mode) {
	switch(mode) {
		case 'r':
		case 'w':
		default:
			cerr << "In function 'openFile': Only accepted modes are 'w' and 'r'.\n"; throw 0;
	}
	return fileName;
}


/**
 * Read command-line arguments
 */
static void parseOptions(int argc, char **argv) {
	int option_index = 0;
	int next_option;
	do {
		next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
		switch (next_option) {
			case ARG_WRAPPER:
				wrapper = optarg;
				break;
			case ARG_USAGE:
			case 'h':
				printUsage(cout);
				throw 0;
				break;
			case 'v': verbose = true; break;
			case ARG_VERSION: showVersion = true; break;
			case ARG_EXCL_AMBIG: exclAllGaps = true; break;
			case ARG_EXTRA: extra = true; break;
			case 'e': refFromEbwt = true; break;
			case 'n': names_only = true; break;
			case 's': summarize_only = true; break;
			case 'a': across = parseInt(-1, "-a/--across arg must be at least 1"); break;
			case 'i': startIndex = parseInt(0, "-i/--start-index must be at least 0"); break;
			case 'j': stopIndex = parseInt(0, "-j/--stop-index must be at least 0"); break;
			case 'c': chrName = optarg;  break;
			case 'f': 
				batchMode = true; 
				resultFileName = optarg; 
				break;
			case 'o': outputFileName = optarg; break;
			case -1: break; /* Done with options. */
			case 0:
				if (long_options[option_index].flag != 0)
					break;
			default:
				printUsage(cerr);
				throw 1;
		}
	} while(next_option != -1);
}

/* This is a little flimsy: Find the header based on the first column: GeneName
 * Return the header line string, with all extra columns appended. 
 * Return empty string if header line wasn't found. */
std::string seekHeaderLine( std::ifstream &resultFile, std::ostream &writeFile ) {
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
				buf = buf.append("\tbracketed sequence");
				buf = buf.append("\tspliced sequence");
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


//getRefSearchData(tokens, chrName, rangeL, rangeR, wideL, wideR, spliceLen);
void getRefSearchData(vector<string>tokens, 
					string &chrName, long &rangeL, 
					long &rangeR, long &wideL, 
					long &wideR, long &spliceLen) 
{
	int tokNo = 0;
	// In this iteration through a line's data sections, grab information but do not 
	// print or edit contents of the tokens
	for(std::vector<std::string>::iterator itr = tokens.begin(); itr != tokens.end(); ++itr ) {
		if( tokNo == CHR_COL ) {
			chrName = std::string(*itr);
		}
		// Since we are just looking at 15 bp on either side of the spliced
		// sequence, I don't think we need to do anything here. -AMC
		// The work can be done by simply +/- BOUNDARY_LEN to rangeL/R
		/*if( tokNo == BRACKET_SEQ_COL ) {

		}*/
		if( tokNo == RANGE_COL ) {
			int pos = itr->find("--");
			rangeL = atol( itr->substr(0,pos).c_str() );
			rangeR = atol( itr->substr(pos+2).c_str() );
			wideL = rangeL - BOUNDARY_LEN;
			wideR = rangeR + BOUNDARY_LEN;
		}
		if( tokNo == SPLICE_LEN_COL ) {
			spliceLen = atol( itr->c_str() );
		}
		tokNo += 1;
	}
	return;
}

void print_fasta_record(ostream& fout,
						const string& defline,
						const string& seq)
{
	fout << ">";
	fout << defline << endl;

	if(across > 0) {
		size_t i = 0;
		while (i + across < seq.length())
		{
			fout << seq.substr(i, across) << endl;
			i += across;
		}
		if (i < seq.length())
			fout << seq.substr(i) << endl;
	} else {
		fout << seq << endl;
	}
}

/**
 * Given output stream, name and length, print a string of Ns with the
 * appropriate number of columns.
 */
void print_alln_ref_sequence(
	ostream& fout,
	const string& name,
	size_t len)
{
	fout << ">" << name << "\n";
	size_t j = 0;
	for(size_t i = 0; i < len; i += across) {
		while(j < len && j < i+across) {
			fout << 'N';
			j++;
		}
		fout << "\n";
	}
}

/**
 * Given output stream, BitPairReference, reference index, name and
 * length, print the whole nucleotide reference with the appropriate
 * number of columns.
 */
/*void print_ref_sequence(
	ostream& fout,
	BitPairReference& ref,
	const string& name,
	size_t refi,
	size_t len)
{
	bool newlines = across > 0;
	int myacross = across > 0 ? across : 60;
	size_t incr = myacross * 1000;
	uint32_t *buf = new uint32_t[(incr + 128)/4];
	fout << ">" << name << "\n";
	for(size_t i = 0; i < len; i += incr) {
		size_t amt = min(incr, len-i);
		assert_leq(amt, incr);
		int off = ref.getStretch(buf, refi, i, amt);
		uint8_t *cb = ((uint8_t*)buf) + off;
		for(size_t j = 0; j < amt; j++) {
			if(newlines && j > 0 && (j % myacross) == 0) fout << "\n";
			assert_range(0, 4, (int)cb[j]);
			fout << "ACGTN"[(int)cb[j]];
		}
		fout << "\n";
	}
	delete buf;
}
*/

/**
 * Given output stream, BitPairReference, reference index, name and
 * length, print the reference from start to stop.
 */
void print_ref_sequence_RSR(
	ostream& fout,
	BitPairReference& ref,
	const string& name,
	size_t refi,
	size_t len)
{
  size_t start, stop;
  if (startIndex > 0) start = startIndex;
  else start = 0;
  
  if (stopIndex > 0) stop = stopIndex;
  else stop = len;

  //  cout << "start = " << start << ", stop = " << stop << endl;
  //  cout << "name = " << name << endl;
  
	bool newlines = across > 0;
	int myacross = across > 0 ? across : 60;
	size_t incr = myacross * 1000;
	uint32_t *buf = new uint32_t[(incr + 128)/4];
	fout << ">" << name << "\n";
	if (stop > len) stop = len;
	for(size_t i = start; i < stop; i += incr) {
		size_t amt = min(incr, stop-i);
		assert_leq(amt, incr);
		int off = ref.getStretch(buf, refi, i, amt);
		uint8_t *cb = ((uint8_t*)buf) + off;
		for(size_t j = 0; j < amt; j++) {
			if(newlines && j > 0 && (j % myacross) == 0) fout << "\n";
			assert_range(0, 4, (int)cb[j]);
			fout << "ACGTN"[(int)cb[j]];
		}
		fout << "\n";
	}
	delete buf;
}

/**
 * Given output stream, BitPairReference, reference index, name and
 * length, return the reference from start to stop in a string.
 */
std::string get_ref_sequence_RSR_string(
	BitPairReference& ref,
	const string& name,
	size_t refi,
	size_t len)
{
	int size = 0;
  size_t start, stop;
	string returnStr = "";
  if (startIndex > 0) start = startIndex;
  else start = 0;
  
  if (stopIndex > 0) stop = stopIndex;
  else stop = len;

	bool newlines = across > 0;
	int myacross = across > 0 ? across : 60;
	size_t incr = myacross * 1000;
	uint32_t *buf = new uint32_t[(incr + 128)/4];
	//fout << ">" << name << "\n";
	if (stop > len) stop = len;
	char cbuf[stop - start + 2];
	for(size_t i = start; i < stop; i += incr) {
		size_t amt = min(incr, stop-i);
		assert_leq(amt, incr);
		int off = ref.getStretch(buf, refi, i, amt);
		//uint8_t i*cb = ((uint8_t*)buf) + off;
		uint8_t *cb = ((uint8_t*)buf) + off;
		for(size_t j = 0; j < amt; j++) {
			//if(newlines && j > 0 && (j % myacross) == 0) fout << "\n";
			assert_range(0, 4, (int)cb[j]);
			//fout << "ACGTN"[(int)cb[j]];
			//returnStr.append("ACGTN"[(int)cb[j]]);
			cbuf[size++] = "ACGTN"[(int)cb[j]];
		}
		//fout << "\n";
	}
	cbuf[size] = '\0';
	returnStr = returnStr.append(cbuf);

	delete buf;
	return returnStr;
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
	long tmp = atol(range.substr(0,endFirstNo).c_str()) - 1; 
	char buf[50];
	snprintf(buf, 50, "%ld", tmp);
	newRange = newRange.append( buf );
	newRange = newRange.append( range.substr(endFirstNo) );
	return newRange;
}

//TODO: adapt this function to work with current format
/*void makeNewStrs( std::string &newBrackSeq, std::string &newSpliceSeq, ile, int spliceRangeLen, int boundaryLen ) {
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
}*/

void print_batch_line(ostream &fout, vector<string> &tokens, string &splSeq, int spliceSeqLen, int boundaryLen) {
	//TODO: makeNewStrs edit
	int tokNo = 0;
	string lBrack, rBrack;

	lBrack = splSeq.substr(0,boundaryLen);
	rBrack = splSeq.substr(boundaryLen + spliceSeqLen);
	splSeq = splSeq.substr(boundaryLen, spliceSeqLen);
	//print the modified current line to the writefile
	for(std::vector<std::string>::iterator itr = tokens.begin(); itr != tokens.end(); ++itr ) {
		if( tokNo == RANGE_COL ) {
			// Edit the range of supp. reads column to show proper inclusion w/ square brackets
			*itr = editRange(*itr);
		}
		fout << *itr << '\t';
		fout.flush();
		tokNo += 1;
	}
	// First, output the bracketed sequence
	string newBrackSeq = lBrack.append("]--[").append(rBrack);
	fout << newBrackSeq << '\t';
	// Then, output the junction site region
	fout << splSeq << endl;
	fout.flush();
	return;
}

/**
 * Create a BitPairReference encapsulating the reference portion of the
 * index at the given basename.  Iterate through the reference
 * sequences, sending each one to print_ref_sequence to print.
 */
void print_ref_sequences(
	ifstream& resultFile,
	ostream& fout,
	bool color,
	const vector<string>& refnames,
	const TIndexOffU* plen,
	const string& adjustedEbwtFileBase)
{
	BitPairReference ref(
		adjustedEbwtFileBase, // input basename
		color,                // true -> expect colorspace reference
		false,                // sanity-check reference
		NULL,                 // infiles
		NULL,                 // originals
		false,                // infiles are sequences
		true,                 // load sequence
		false,                // memory-map
		false,                // use shared memory
		false,                // sweep mm-mapped ref
		verbose,              // be talkative
		verbose);             // be talkative at startup






#ifdef ACCOUNT_FOR_ALL_GAP_REFS
	if( batchMode ) {
		long wideL, wideR, rangeL, rangeR, spliceLen;
		string buf, bracketL, bracketR, splSeq;
		vector<string> tokens;

		getline(resultFile, buf); 
		while( !resultFile.eof() ) {
			tokens = splitString(buf, '\t');
			getRefSearchData(tokens, chrName, rangeL, rangeR, wideL, wideR, spliceLen); 
			startIndex = wideL;
			stopIndex = wideR;

			for(size_t i = 0; i < ref.numNonGapRefs(); i++) {
				if( chrName == refnames[i] ) {
					splSeq = get_ref_sequence_RSR_string(
						ref,
						chrName,
						i,
						plen[i] + (color ? 1 : 0));
				}
			}
			//splSeq contains the spliced sequence + 15 extra nucl. on each side
			print_batch_line(fout, tokens, splSeq, rangeR-rangeL, BOUNDARY_LEN );  
			getline(resultFile, buf); 
		}
	}
	else {	//not batch mode
		for(size_t i = 0; i < ref.numRefs(); i++) {
			if(ref.isAllGaps(i) && !exclAllGaps) {
				if (chrName == "" || chrName == refnames[i])
				print_alln_ref_sequence_RSR(
					fout,
					refnames[i],
					ref.len(i));
			} else {
				if (chrName == "" || chrName == refnames[i])
				print_ref_sequence_RSR(
					fout,
					ref,
					refnames[i],
					ref.shrinkIdx(i),
					ref.len(i));
			}
		}
	}
#else
	assert_eq(refnames.size(), ref.numNonGapRefs());
	if( batchMode ) {
		long wideL, wideR, rangeL, rangeR, spliceLen;
		string buf, bracketL, bracketR, splSeq;
		vector<string> tokens;

		getline(resultFile, buf); 
		while( !resultFile.eof() ) {
			tokens = splitString(buf, '\t');
			getRefSearchData(tokens, chrName, rangeL, rangeR, wideL, wideR, spliceLen); 
			startIndex = wideL;
			stopIndex = wideR;

			for(size_t i = 0; i < ref.numNonGapRefs(); i++) {
				if( chrName == refnames[i] ) {
					splSeq = get_ref_sequence_RSR_string(
						ref,
						chrName,
						i,
						plen[i] + (color ? 1 : 0));
				}
			}
			//splSeq contains the spliced sequence + BOUNDARY_LEN  extra nucl. on each side
			print_batch_line(fout, tokens, splSeq, rangeR-rangeL, BOUNDARY_LEN );  
			getline(resultFile, buf); 
		}
	}
	else { //not batch mode 
		for(size_t i = 0; i < ref.numNonGapRefs(); i++) {
			if (chrName == "" || chrName == refnames[i]) {
				print_ref_sequence_RSR(
					fout,
					ref,
					refnames[i],
					i,
					plen[i] + (color ? 1 : 0));
			}
		}
	}
#endif
}

/**
 * Given an index, reconstruct the reference by LF mapping through the
 * entire thing.
 */
template<typename TStr>
void print_index_sequences(
	ostream& fout,
	Ebwt<TStr>& ebwt,
	const BitPairReference& refs)
{
	vector<string>* refnames = &(ebwt.refnames());

	TStr cat_ref;
	ebwt.restore(cat_ref);

	TIndexOffU curr_ref = OFF_MASK;
	string curr_ref_seq = "";
	TIndexOffU curr_ref_len = OFF_MASK;
	uint32_t last_text_off = 0;
	size_t orig_len = seqan::length(cat_ref);
	TIndexOffU tlen = OFF_MASK;
	bool first = true;
	for(size_t i = 0; i < orig_len; i++) {
		TIndexOffU tidx = OFF_MASK;
		TIndexOffU textoff = OFF_MASK;
		tlen = OFF_MASK;

		ebwt.joinedToTextOff(1 /* qlen */, (TIndexOffU)i, tidx, textoff, tlen);

		if (tidx != OFF_MASK && textoff < tlen)
		{
			if (curr_ref != tidx)
			{
				if (curr_ref != OFF_MASK)
				{
					// Add trailing gaps, if any exist
					if(curr_ref_seq.length() < curr_ref_len) {
						curr_ref_seq += string(curr_ref_len - curr_ref_seq.length(), 'N');
					}
					print_fasta_record(fout, (*refnames)[curr_ref], curr_ref_seq);
				}
				curr_ref = tidx;
				curr_ref_seq = "";
				curr_ref_len = tlen;
				last_text_off = 0;
				first = true;
			}

			TIndexOffU textoff_adj = textoff;
			if(first && textoff > 0) textoff_adj++;
			if (textoff_adj - last_text_off > 1)
				curr_ref_seq += string(textoff_adj - last_text_off - 1, 'N');

			curr_ref_seq.push_back(getValue(cat_ref,i));
			last_text_off = textoff;
			first = false;
		}
	}
	if (curr_ref < refnames->size())
	{
		// Add trailing gaps, if any exist
		if(curr_ref_seq.length() < curr_ref_len) {
			curr_ref_seq += string(curr_ref_len - curr_ref_seq.length(), 'N');
		}
		print_fasta_record(fout, (*refnames)[curr_ref], curr_ref_seq);
	}

}

static char *argv0 = NULL;

void print_index_sequence_names(const string& fname, ostream& fout)
{
	vector<string> p_refnames;
	readEbwtRefnames(fname, p_refnames);
	for(size_t i = 0; i < p_refnames.size(); i++) {
		cout << p_refnames[i] << endl;
	}
}

typedef Ebwt<String<Dna, Packed<Alloc<> > > > TPackedEbwt;

/**
 * Print a short summary of what's in the index and its flags.
 */
void print_index_summary(
	const string& fname,
	ostream& fout,
	const BitPairReference& refs)
{
	int32_t flags = readFlags(fname);
	int32_t flagsr = readFlags(fname + ".rev");
	bool color = readEbwtColor(fname);
	bool entireReverse = readEntireReverse(fname + ".rev");
	TPackedEbwt ebwt(
		fname,
		color,                // index is colorspace
		-1,                   // don't require entire reverse
		true,                 // index is for the forward direction
		-1,                   // offrate (-1 = index default)
		-1,
		false,                // use memory-mapped IO
		false,                // use shared memory
		false,                // sweep memory-mapped memory
		true,                 // load names?
		//false,                // load SA sample?
		NULL,                 // no reference map
		verbose,              // be talkative?
		verbose,              // be talkative at startup?
		false,                // pass up memory exceptions?
		false);               // sanity check?
	vector<string> p_refnames;
	readEbwtRefnames(fname, p_refnames);
	if(extra) {
		cout << "Flags" << '\t' << (-flags) << endl;
		cout << "Reverse flags" << '\t' << (-flagsr) << endl;
	}
	cout << "Colorspace" << '\t' << (color ? "1" : "0") << endl;
	if(extra) {
		cout << "Concat then reverse" << '\t' << (entireReverse ? "1" : "0") << endl;
		cout << "Reverse then concat" << '\t' << (entireReverse ? "0" : "1") << endl;
		cout << "nPat" << '\t' << ebwt.nPat() << endl;
		cout << "refnames.size()" << '\t' << p_refnames.size() << endl;
		cout << "refs.numRefs()" << '\t' << refs.numRefs() << endl;
		cout << "refs.numNonGapRefs()" << '\t' << refs.numNonGapRefs() << endl;
	}
	cout << "SA-Sample" << "\t1 in " << (1 << ebwt.eh().offRate()) << endl;
	cout << "FTab-Chars" << '\t' << ebwt.eh().ftabChars() << endl;
	for(size_t i = 0; i < ebwt.nPat(); i++) {
		cout << "Sequence-" << (i+1)
		     << '\t' << p_refnames[refs.expandIdx((uint32_t)i)]
		     << '\t' << (ebwt.plen()[i] + (color ? 1 : 0))
		     << endl;
	}
	if(extra) {
		cout << "RefRecords:\n";
		for(size_t i = 0; i < refs.refRecords().size(); i++) {
			RefRecord r = refs.refRecords()[i];
			cout << r.first << "\t(" << r.off << ", " << r.len << ")" << endl;
		}
	}
}

static void driver(
	const string& ebwtFileBase,
	const string& query)
{
	// Adjust
	string adjustedEbwtFileBase = adjustEbwtBase(argv0, ebwtFileBase, verbose);
	if (names_only) {
		print_index_sequence_names(adjustedEbwtFileBase, cout);
		return;
	}
	bool color = readEbwtColor(adjustedEbwtFileBase);
	BitPairReference refs(
		adjustedEbwtFileBase,
		color,
		false,
		NULL,
		NULL,
		false,
		false, // don't load sequence (yet)
		false,
		false,
		false, // mmSweep
		verbose,
		verbose);
	if(summarize_only) {
		print_index_summary(adjustedEbwtFileBase, cout, refs);
	} else {
		// Initialize Ebwt object
		TPackedEbwt ebwt(
			adjustedEbwtFileBase,
			color,                // index is colorspace
			-1,                   // don't care about entire-reverse
			true,                 // index is for the forward direction
			-1,                   // offrate (-1 = index default)
			-1,
			false,                // use memory-mapped IO
			false,                // use shared memory
			false,                // sweep memory-mapped memory
			true,                 // load names?
			//true,                 // load SA sample?
			NULL,                 // no reference map
			verbose,              // be talkative?
			verbose,              // be talkative at startup?
			false,                // pass up memory exceptions?
			false);               // sanity check?
		// Load whole index into memory
		if(refFromEbwt) {
			ebwt.loadIntoMemory(-1, -1, true, false);
			print_index_sequences(cout, ebwt, refs);
		} else {
			vector<string> refnames;
			readEbwtRefnames(adjustedEbwtFileBase, refnames);
			// You need an output stream regardless
			ifstream resultFile;
			ofstream writeFile;
			// The following 2 if/[else] blocks handle file I/O
			if( !resultFileName.empty() ) {
				resultFile.open(resultFileName.c_str(), std::ifstream::in);
				if( !resultFile.good() ) {
					cerr << "Error, could not open result file '" << resultFileName << "'.";
					cerr << endl << "--->Does it exist?" << endl;
					throw 1;
				}
			}
			// Output to standard output unless an output file is specified
			ostream &outFile = (outputFileName.empty() ? cout : writeFile );
			if( !outputFileName.empty() ) {
				// Choose between default filename....
				if( outputFileName == string("default") ) {
					if( resultFileName.empty() ) {
						cerr << "To use default output option, include result file." << endl;
						printUsage(cerr);
						throw 1;
					}
					writeFile.open(resultFileName.append(".splSeq").c_str(), std::ofstream::out); 
				} 
				// ....or provided filename
				else {
					writeFile.open(outputFileName.c_str(), std::ofstream::out); 
				}
				if( !writeFile.good() ) {
					cerr << "Error opening output file. '" << outputFileName << "'." << endl;
					throw 1;
				}
			}
			if( batchMode ) {
				seekHeaderLine(resultFile, outFile);
			}
			print_ref_sequences(
				resultFile,
				outFile,
				readEbwtColor(ebwtFileBase),
				refnames,
				ebwt.plen(),
				adjustedEbwtFileBase);
		}
		// Evict any loaded indexes from memory
		if(ebwt.isInMemory()) {
			ebwt.evictFromMemory();
		}
	}
}

/**
 * main function.  Parses command-line arguments.
 */
int main(int argc, char **argv) {
	try {
		string ebwtFile;  // read serialized Ebwt from this file
		string query;   // read query string(s) from this file
		vector<string> queries;
		string outfile; // write query results to this file
		argv0 = argv[0];
		parseOptions(argc, argv);
		if(showVersion) {
			cout << argv0 << " version " << BOWTIE_VERSION << endl;
			if(sizeof(void*) == 4) {
				cout << "32-bit" << endl;
			} else if(sizeof(void*) == 8) {
				cout << "64-bit" << endl;
			} else {
				cout << "Neither 32- nor 64-bit: sizeof(void*) = " << sizeof(void*) << endl;
			}
			cout << "Built on " << BUILD_HOST << endl;
			cout << BUILD_TIME << endl;
			cout << "Compiler: " << COMPILER_VERSION << endl;
			cout << "Options: " << COMPILER_OPTIONS << endl;
			cout << "Sizeof {int, long, long long, void*, size_t, off_t}: {"
				 << sizeof(int)
				 << ", " << sizeof(long) << ", " << sizeof(long long)
				 << ", " << sizeof(void *) << ", " << sizeof(size_t)
				 << ", " << sizeof(off_t) << "}" << endl;
			return 0;
		}

		// Get input filename
		if(optind >= argc) {
			cerr << "No index name given!" << endl;
			printUsage(cerr);
			return 1;
		}
		ebwtFile = argv[optind++];

		// Optionally summarize
		if(verbose) {
			cout << "Input ebwt file: \"" << ebwtFile << "\"" << endl;
			cout << "Output file: \"" << outfile << "\"" << endl;
			cout << "Local endianness: " << (currentlyBigEndian()? "big":"little") << endl;
#ifdef NDEBUG
			cout << "Assertions: disabled" << endl;
#else
			cout << "Assertions: enabled" << endl;
#endif
		}
		driver(ebwtFile, query);
		return 0;
	} catch(std::exception& e) {
		cerr << "Command: ";
		for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
		cerr << endl;
		return 1;
	} catch(int e) {
		if(e != 0) {
			cerr << "Command: ";
			for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
			cerr << endl;
		}
		return e;
	}
}
