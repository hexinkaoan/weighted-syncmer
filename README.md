weighted-syncmer
========================================================================

We introduce weighted syncmer sampling, which integrates weighted minimizer sampling and syncmer sampling to enhance the sensitivity and accuracy of long-read mapping. We modified two state-of-the-art long-read mappers, Minimap2 and Winnowmap, by substituting the sketching sampling methods with weighted syncmer sampling. We assessed their sensitivity and accuracy using simulated and real datasets. The experimental results indicate that weighted syncmer sampling significantly improves the sensitivity and accuracy of long-read mapping. Detailed evaluations are available from the paper [Weighted syncmer sampling improves long-read mapping](https://ieeexplore.ieee.org/abstract/document/10821776)

## Compile

Clone source code from master branch.
  ```sh
	git clone https://github.com/hexinkaoan/weighted-syncmer.git
  ```
[ws-winnowmap](https://github.com/hexinkaoan/weighted-syncmer/tree/main/ws-winnowmap) compilation requires C++ compiler with c++11 and openmp, which are available by default in GCC >= 4.8. Expect `winnowmap` and `meryl` executables in `bin` folder.

  ```sh
	cd ws-winnowmap
	make -j8
  ```
[ws-minimap2](https://github.com/hexinkaoan/weighted-syncmer/tree/main/ws-minimap2) need to have a C compiler, GNU make and zlib development files installed. Then type `make` in the source code directory to compile. 

  ```sh
	cd ws-minimap2
	make
  ```

## Usage

For either mapping long reads or computing whole-genome alignments, weighted-syncmer requires pre-computing high frequency k-mers (e.g., top 0.02% most frequent) in a reference. weighted-syncmer uses [meryl](https://github.com/marbl/meryl) k-mer counting tool for this purpose.  

*  ws-winnowmap
```sh
meryl count k=19 output merylDB ref.fa
meryl print greater-than distinct=0.9998 merylDB > repetitive_k19.txt
	  
./bin/winnowmap -k 19 --syncs 15 --synct 3  -W repetitive_k19.txt -a ref.fa query.fq > alignment.sam
```

*  ws-minimap2

```sh
meryl count k=15 output merylDB ref.fa
meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
	  
./minimap2 -k 19 --syncs 15 --synct 3  -W repetitive_k15.txt -a ref.fa query.fq > out.sam
```

### Algorithm overview

The weighted syncmer sampling method improves the weighting function part in the syncmer sampling method, similar to the weighted minimizer method. This method re-assigns weights based on whether the s-mer is a highly repetitive sequence after assigning weights by the weight hash function.

<div align="center">
  <img src="https://i.postimg.cc/Xqv4qyMw/fig1a.png" width="400px"><br>
  (a) syncmer sampling
</div>


<div align="center">
  <img src="https://i.postimg.cc/7h78f1sJ/fig1b.png" width="400px"><br>
  (b) weighted-syncmer sampling
</div>

## Citations

- **Chirag Jain, Arang Rhie, Haowen Zhang, Chaudia Chu, Brian Walenz, Sergey Koren and Adam Phillippy**. "[Weighted minimizer sampling improves long read mapping](https://doi.org/10.1093/bioinformatics/btaa435)". *Bioinformatics (ISMB proceedings)*, 2020.
- **Li, H**. "[Minimap2: pairwise alignment for nucleotide sequences.](https://doi.org/10.1093/bioinformatics/bty191)". *Bioinformatics*, 34:3094-3100. doi:10.1093/bioinformatics/bty191
- **Shaw, Jim and Yu, Yun William**. "[Theory of local k-mer selection with applications to long-read alignment.](https://doi.org/10.1093/bioinformatics/btab790)". *Bioinformatics(2022)*
- **Wei Quan; Jinjun Kang; Yanfei Deng; Zhuang Liu; Xiao Zhu; Guangri Quan**. "[[Weighted syncmer sampling improves long-read mapping](https://ieeexplore.ieee.org/abstract/document/10821776)]".
