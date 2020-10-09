change log:
—————————————————————
version january 11, 2017
	+ bug-fix in overlay-edits.awk where id and color fields were switched
—————————————————————
version december 22, 2016
	+ made labels to be applied variables (default is label1=":::fragment_", label2=":::debris)
—————————————————————
version december 14, 2016
	+ renamed run-misassembly-detector.sh to run-mismatch-detector.sh
	+ updated USAGE to highlight thinning
	+ substituted bc calc for threshold with awk handling for consistency
—————————————————————
version december 13, 2016
	+ in overlay-edits.sh added separation between intra-scaffold and inter-scaffold misassemblies.
—————————————————————
version december 12, 2016
	+ added edge-thinning as one-liner in run-misassembly-detector.sh.
—————————————————————
version december 3, 2016
	+ first version
—————————————————————
This is a group of scripts to highlight and excise misassemblies.