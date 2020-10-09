change log:
version december 14, 2016
	+ substituted juicebox_tools.jar for a version with working -r and -n flags
	+ made appropriate modifications in juicebox_tools.sh
	+ added explicit building of the options string based on zoom
	+ added unprompted option -n for no normalization
	+ added unprompted -c flag for clean-up after done to delete the lifted mnd file
—————————————————————
version december 7, 2016.
	+ renamed script
	+ substitute explicit paths to dependencies to check in same folder as wrapper script
	+ remove option to pass genome id
	+ change argument order
	+ remove option to pass stats and graph files
	+ added unscaled outputs to superscaffold track output
	+ fixed superscaffold annotation 1-offset
—————————————————————
version december 1, 2016
	+ substituted true fragments for mock fragments. TODO: do proper fragment number remapping to achieve inner peace
—————————————————————
version july 18, 2016
	+ added GNU Parallel wrapper as an option
—————————————————————
earlier versions
—————————————————————
This is a script to visualize draft assemblies based on Juicebox pre function from Juicebox Command Line Tools.