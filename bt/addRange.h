
#include <vector>
#include <string>
#include <fstream>

const std::string TMP_FILENAME = ".splSeq";
const std::string HEADER_ID = "genename";
const int BOUNDARY_LEN = 15;


const int GENE_COL = 0;
const int CHR_COL = 1;
const int SUPPORT_READS_COL = 2;
const int SUPPORT_HALF_COL = 3;
const int SUPPORT_TOTAL_COL = 4;
const int SPLICE_LEN_COL = 5;
const int RANGE_COL = 6;
const int NOVEL_COL = 7;
const int BRACKET_SEQ_COL = 8;
const int SPLICED_SEQ_COL = 9;
 

void usage();
void bracketRange(std::string &range);
void modifyResults(std::ifstream& resultFile, std::ofstream& writeFile, std::string genFileName);
void extendRange(std::string token, std::string &bracketL, std::string &bracketR, long &moveLeft, long &moveRight);
//int getColumnNo( std::string header, std::string column );
void makeNewStrs( std::string &newBrackSeq, std::string &newSpliceSeq, FILE *spliceFile, int spliceLen, int boundaryLen );
/*void makeNewStrs_saveEdges( std::string &newBrackSeq, std::string &newSpliceSeq, FILE *spliceFile, int bp, long wideL,
                  long rangeL, long rangeR, long wideR, std::string &bracketL, std::string &bracketR );*/
std::string seekHeaderLine( std::ifstream& resultFile, std::ostream& writeFile );
std::string editRange(std::string& range);
std::vector<std::string> splitString(std::string str, const char delimit);

