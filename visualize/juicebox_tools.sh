#!/bin/sh
set -e
java -Xms49152m -Xmx49152m -jar `dirname $0`/juicebox_tools.jar $*
