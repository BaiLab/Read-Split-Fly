#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <libgen.h>

const unsigned int K = 1024;
int MIN_SPLIT_LENGTH;
const char *SUFFIX = ".split";

void parse(FILE *in, FILE *out) {
    //char id[K],seq[K],str[K],score[K], tmp[K];
    char *id,*seq,*str,*score,tmp;
    size_t n;
    long i = 0u, k = 0u, l = 0u;
    int j;
    int left_cut, right_cut; 
    id = seq = str = score = 0;
    while(!feof(in)) {
        left_cut = MIN_SPLIT_LENGTH, right_cut = left_cut;

        if (getline(&id, &n, in) == -1) break;
        if (id[strlen(id)-1]== '\n') id[strlen(id)-1] = 0;
        if (getline(&seq, &n, in) == -1) break;
        if (seq[strlen(seq)-1] == '\n') seq[strlen(seq)-1] = 0;
        const int readLength = strlen(seq);

        right_cut = readLength - MIN_SPLIT_LENGTH; // Jeff, to handle different length sequences

        if (getline(&str, &n, in) == -1) break;
        if (str[strlen(str)-1] == '\n') str[strlen(str)-1] = 0;

        if (getline(&score, &n, in) == -1) break;
        if (score[strlen(score)-1] == '\n') score[strlen(score)-1] = 0;

        k+=4;
        l = 0u;
        //must have all 4 to contnue...
        while(left_cut <= right_cut) {
            ++l;

            // note: also now printing out the total length of the sequence, because that can vary from read to read for some datasets.
    
            //left
            fprintf(out, "%s-L-%d-%d\n",id,left_cut,readLength);
            fprintf(out, "%.*s\n",left_cut,seq);
            fprintf(out, "%s-L-%d-%d\n",str,left_cut,readLength);
            fprintf(out, "%.*s\n",left_cut,score);
    
            //right
            fprintf(out, "%s-R-%d-%d\n",id,right_cut,readLength);
            fprintf(out, "%.*s\n",right_cut,seq+left_cut);
            fprintf(out, "%s-R-%d-%d\n",str,right_cut,readLength);
            fprintf(out, "%.*s\n",right_cut,score+left_cut);

            //REVERSED!
            //left
            fprintf(out, "%s-L-%d-%d\n",id,right_cut,readLength);
            fprintf(out, "%.*s\n",right_cut,seq);
            fprintf(out, "%s-L-%d-%d\n",str,right_cut,readLength);
            fprintf(out, "%.*s\n",right_cut, score);
 
            //right
            fprintf(out, "%s-R-%d-%d\n",id,left_cut,readLength);
            fprintf(out, "%.*s\n",left_cut,seq+right_cut);
            fprintf(out, "%s-R-%d-%d\n",str,left_cut,readLength);
            fprintf(out, "%.*s\n",left_cut,score+right_cut);

            //next
            ++left_cut; --right_cut;
            i+=16;
        }
        free(id); id = 0;
        free(seq); seq = 0;
        free(str); str = 0;
        free(score); score = 0;
    }
    fprintf(stderr,"in: %lu, out: %lu, passes each line: %lu, total passes: %lu\n", k, i, l, l*k);
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr,"Usage; %s <file> <min size>\n",argv[0]);
        return 1;
    }
    clock();
    clock_t after, before;
    before = clock();
    FILE *fr, *fw;
    char *newfn, *name;
    int sp = 0, sn, ss = strlen(SUFFIX);
    char flag = 0;
    if (argc >= 4) {
        sp = strlen(argv[3]);
        if ((flag = ('/' == argv[3][sp-1]))) ++sp;
        name = basename(argv[1]);
    } else {
        name = argv[1];
    }
    sn = strlen(name);
    newfn = (char *)malloc(sizeof(char) * (sp+1+sn+ss+1)); // bug fix 6 jan 2016, added +1 to include space for terminating NULL character.  and 28 may 2016, another +1 for "/"
    memset(newfn,0,sizeof(char)*(sp+1+sn+ss+1)); // bug fix 6 jan 2016, added +1 to include space for terminating NULL character
    if (sp) strcpy(newfn,argv[3]);
    if (sp && !flag) strcat(newfn,"/");
    newfn = strcat(strcat(newfn, name),SUFFIX);

    MIN_SPLIT_LENGTH = atoi(argv[2]);

    if (!(fr = fopen(argv[1], "r"))) {
        fprintf(stderr,"Unable to open %s for reading.  ", argv[1]);
        perror("");
        return 1;
    }

    if (!(fw = fopen(newfn, "w"))) {
        fprintf(stderr,"Unable to open %s for reading.  ", argv[1]);
        perror("");
        fclose(fr);
        return 1;
    }
    
    parse(fr,fw);
    fclose(fr);
    fclose(fw);
    after = clock();
    fprintf(stderr,"time elapsed: %lf seconds.\n",(after-before)*(1.0/CLOCKS_PER_SEC));
    return 0;
}
