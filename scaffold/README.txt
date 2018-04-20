change log:
—————————————————————
version may 18, 2017.
	+ some safeguards when merging scores in tiger and liger
	+ including proper handling of remaining div by zero scenario
	+ sort -n substituted for sort -g, parallel sort parameters changed in liger and tiger
—————————————————————
version dec 28, 2016.
	+ added a merge-scores line for no parallel option as the normalization no longer happens inside the scrape script
	+ updated TIGer code
—————————————————————
version dec 12, 2016.
	+ add minimal contig/scaffold size option
	+ default orientation flag set to false
	+ use GNU Parallel is made optional with a -p flag. Check added to see if in path
	+ parallel sort for confidence file
	+ uncomment cleanup of helper files
	+ separate from visualization
	+ remove option to pass a list of mnd files: rare scenario not well suited for parallelization and editing
	+ changed output name to .asm instead of .final.asm
	+ added a check at loop entry to break if input is a single scaffold
	
	+ TODO: safe stop wo div by zero
	+ TODO: safe misassembly annotation (now name~/debris$/)
	+ TODO: maybe get rid of cprops in favor of chrom.sizes styled file. Not more more to it.
—————————————————————
version july 29, 2016.
change log:
	+ added score file scraping parallelization. Requires GNU parallel.
	+ fixed a read position determination issue and position out of bound problem.
	+ joined orientation LIGer and regular LIGer. Use -o true/false option to call one of the two. Orientation liger is the current default.
version may 3, 2015.
change log:
	+ a check for algorithm stalling is added as a separate subscript and as part of liger. This identifies problematic contigs and sets them aside starting from the smallest. The required script drop-smallest-dubious-elemt.awk is added and the appropriate changes in the main wrapper script are made
	+ criteria for algorithm finishing has been changed to incorporate the anti-stalling portion. These will likely have to be revisited as not all of the scenarios are properly handled
—————————————————————
version april 23, 2015.
change log:
	+ input changed: instead of the fasta file cprops file now expected for logistic ease
	+ generate-cprops-file script no longer necessary
	+ option -s associated with contig size filtering are removed - assumed this is done upstream during the generation of the cprops file
	+ thresholding in accept-links script changed such that everything >1 is accepted
	
—————————————————————
version march 20, 2015.

This is a stand-alone version of the iterative contig merging algorithm that is based on contig-wide contact data (as opposed to hash-tag contact data in the original TIGer). This approach is potentially more robust to misassemblies in the original contig set. The potential drawback is overcorrection of contact data for larger contigs and sensitivity to coverage biases. Proper comparison with TIGer to resolve the relative benefits and problems is yet to be performed.