#!/bin/sh
set -e
java -Djava.library.path=/aidenlab/natives/ -Xms49152m -Xmx49152m -jar `dirname $0`/juicebox_tools.jar $*
