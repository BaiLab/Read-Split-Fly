all: sp4 sfc srr sbc compare blast_dir bt_dir

sp4: 
	g++ -O4 -o sp4 src/splitPairs.cpp -std=c++11

sfc: 
	gcc -O4 -o sfc src/split_columns.c

srr: 
	gcc -g -o srr src/split_read_rsw.c

sbc: 
	gcc -O4 -o sbc src/split_on_chrom.c

compare:
	g++ -O4 -o compare src/compare.cpp -std=c++11

blast_dir: blast/makeFASTA
bt_dir: bt/bowtie-inspect-l-RSR

bt/bowtie-inspect-l-RSR: 
	$(MAKE) -C bt

blast/makeFASTA:
	$(MAKE) -C blast

.PHONY: clean-small
clean-small:
	rm -f sp4 sfc srr sbc compare
	$(MAKE) -C blast clean
  
.PHONY: clean
clean: clean-small
	$(MAKE) -C bt clean
  
