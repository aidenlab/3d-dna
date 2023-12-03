# 3D de novo assembly (3D-DNA) pipeline

This version of the pipeline (180419) updates the merge section that aims to address errors of undercollapsed heterozygosity.

For a detailed description of the pipeline and how it integrates with other tools designed by the Aiden Lab see [Genome Assembly Cookbook](http://aidenlab.org/assembly/manual_180322.pdf) on http://aidenlab.org/assembly.

For the original version of the pipeline and to reproduce the Hs2-HiC and the AaegL4 genomes reported in [(Dudchenko et al., *Science*, 2017)](http://science.sciencemag.org/content/356/6333/92) see the [original commit](https://github.com/theaidenlab/3d-dna/tree/745779bdf64db6e55bddb70c24e9b58825938c33).

For the detailed description of the merge section see https://github.com/theaidenlab/AGWG-merge.

Feel free to post your questions and comments at:
[http://www.aidenlab.org/forum.html]

### Citations and licensing
If you use this code or the resulting assemblies, please cite the following paper:

Dudchenko, O., Batra, S.S., Omer, A.D., Nyquist, S.K., Hoeger, M., Durand, N.C., Shamim, M.S., Machol, I., Lander, E.S., Aiden, A.P., et al. (2017). *De novo assembly of the *Aedes aegypti* genome using Hi-C yields chromosome-length scaffolds.* Science. Apr 7; 356(6333):92-95. doi: https://doi.org/10.1126/science.aal3327. Epub 2017 Mar 23.

#### This software is distributed under The MIT License (MIT).

### Overview of the pipeline
`Software version: 180922`

An overview of the detailed workflow of the 3D-DNA pipeline is schematically given in Fig. S1 in [(Dudchenko et al., *Science*, 2017)](http://science.sciencemag.org/content/356/6333/92).

We begin with a series of iterative steps whose goal is to eliminate misjoins in the input scaffolds. Each step begins with a scaffold pool (initially, this pool is the set of input scaffolds themselves). The scaffolding algorithm is used to order and orient these scaffolds. Next, the misjoin correction algorithm is applied to detect errors in the scaffold pool, thus creating an edited scaffold pool. Finally, the edited scaffold pool is used as an input for the next iteration of the misjoin correction algorithm. The ultimate effect of these iterations is to reliably detect misjoins in the input scaffolds without removing correctly assembled sequence. After this process is complete, the scaffolding algorithm is applied to the revised input scaffolds, and the output – a single “megascaffold” which concatenates all the chromosomes – is retained for post-processing.

This post-processing includes four steps: (i) a polishing algorithm which attempts to fix errors when cumulative 3D signal 'wins' over the diagonal; (ii) a chromosome splitting algorithm, which is used to extract the chromosome-length scaffolds from the megascaffold; (iii) a sealing algorithm, which detects false positives in the misjoin correction process, and restores the erroneously removed sequence from the original scaffold; and (iv) a merge algorithm, which corrects misassembly errors due to undercollapsed heterozygosity in the input scaffolds. Step (iv) are omitted for Hs2-HiC, which is not in the Rabl configuration and lacks substantial undercollapsed heterozygosity.

#### Prerequisites
- `LastZ (version 1.03.73 released 20150708)` – for diploid mode only
- `Java version >=1.7`
- `Bash >=4`
- `GNU Awk >=4.0.2`
- `GNU coreutils sort >=8.11`
- `Python >=2.7` -  for chromosome number-aware splitter module only 
- `scipy numpy matplotlib` - for chromosome number-aware splitter module only
- `BioPython` - for edit fasta
- `seqkit` - for wrapping sequence

#### Recommended
- `GNU Parallel >=20150322` – highly recommended to increase performance

#### Installation
There is no need for installation.

#### Overview
The pipeline consists of one bash wrapper script [run-asm-pipeline.sh](https://github.com/theaidenlab/3d-dna/blob/master/run-asm-pipeline.sh) that calls individual modules to assemble a genome. Run the following to list expected arguments: 

```sh
./run-asm-pipeline.sh –h
```
and run
```sh
./run-asm-pipeline.sh –-help
```
to see a full list of available options.

Another script [run-asm-pipeline-post-review.sh](https://github.com/theaidenlab/3d-dna/blob/master/run-asm-pipeline-post-review.sh) is designed to execute the later stages of the pipeline (finalize or seal and finalize) when picking up after Juicebox Assembly Tools (JBAT) review. See [Genome Assembly Cookbook](http://aidenlab.org/assembly/manual_180322.pdf) and also [this preprint](https://www.biorxiv.org/content/early/2018/01/28/254797) for details about JBAT.

Run the following for the description of arguments and options: 
```sh
./run-asm-pipeline-post-review.sh –h|--help
```

### Individual modules
The wrapper script calls individual modules of the pipeline. Code related to any particular module is organized into individual folders. The modules can also be run as separate scripts. The list of individual modules with their main wrapper scripts is given below:

a) `	scaffold (./run-liger-scaffolder.sh)`	– scaffold the draft assembly;

b) `	visualize (./run-asm-visualizer.sh)`	– make Juicebox-compatible contact map;

c) `	edit (./run-misassembly-detector.sh)`	– annotate regions consistent with a misassembly;

d) `	polish (./run-asm-polisher.sh)`	– polish assembly;

e) `	split (./run-asm-splitter.sh)`	– split scaffolder output into raw chromosomal scaffolds. Note that though the current version of the pipeline does not expect to know the number of chromosomes expected a chromosome-number-aware version of the splitter is availabe as a standalone module;

f) `	seal (./seal-asm.sh)`	– reintroduce some fragments back into the assembly;

g) `	merge (./run-asm-merger.sh)`	– merge fragments of the assembly that represent alternative haplotypes (only in diploid mode);

h) `	finalize (./finalize-output.sh)`	– generates an output fasta;

Additional modules include:

i) `	utils`	– several core scripts that are used across modules;

j) `	lift`	– several core scripts to liftover coordinates from the draft to assembly and vice versa;

k) `	supp`	– several additional scripts to match available map data and generate initial conditions;

l) `	data`	– mapping data tables used for validation of the 2017 assemblies (AaegL4, CpipJ3 and Hs2-HiC).


### Output files
The pipeline generates a number of files. The types of files are listed below. The main outputs are the fasta files annotated as “FINAL” which contains the output scaffolds.

a) 	.fasta files
- 	“FINAL” – chromosome-length scaffolds;
- 	“final” – input with all the misjoin correction introduced;

b) 	.hic files
-	“FINAL“ - after the addition of gaps to the chromosome-length assembly (built on request with --build-gapped-map option);
-   	“final“ - after sealing stage;
- 	“polished” – after polishing stage;
- 	“resolved” – after editing and scaffolding;
- 	[0123…] – correspond to the assembly at individual editing iterations;

c) 	.scaffold_track.txt & .superscaf_track.txt
- 	Scaffold boundary files (Juicebox 2D annotation format) - still supported but now relatively obsolete that .assembly files can be loaded into Juicebox via Juicebox Assembly Tools;

d) 	.bed & .wig files
- 	Tracks illustrating putative misjoins;

e) 	.assembly (supersedes .cprops and .asm files)
- 	Custom file format that tracks modifications to the input contigs at various stages in the assembly. Together with matching .hic files input to Juicebox Assembly Tools. The files are available for all stages of the assembly including:
-	“FINAL“ - after the addition of gaps to the chromosome-length assembly;
-   	“final“ - after sealing stage;
- 	“polished” – after polishing stage;
- 	“resolved” – after editing and scaffolding;
- 	[0123…] – correspond to the assembly at individual editing iterations;


f) 	supplementary files:
- 	edits.for.step.\*.txt; mismatches.at.step.\*.txt; suspect_2D.at.step.\*.txt - list of problematic regions (Juicebox 2D annotation format);
- 	alignments.txt (for diploid mode only) – pairwise alignment data for alternative haplotype candidates, parsed from LASTZ output.


[Supporting Online Materials]: <http://science.sciencemag.org/content/suppl/2017/03/22/science.aal3327.DC1?_ga=1.9816115.760837492.1490574064>
[Dudchenko et al., De novo assembly of the Aedes aegypti genome using Hi-C yields chromosome-length scaffolds. Science, 2017.]: <http://science.sciencemag.org/content/early/2017/03/22/science.aal3327.full>
[Juicer (Durand & Shamim et al., Cell Systems, 2016)]: <http://www.cell.com/cell-systems/abstract/S2405-4712(16)30219-8>
[Juicebox (Durand & Robinson et al., Cell Systems, 2016)]: <http://www.cell.com/cell-systems/abstract/S2405-4712(15)00054-X>
[http://www.aidenlab.org/forum.html]: <http://www.aidenlab.org/forum.html>
[GSE95797]: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95797>
