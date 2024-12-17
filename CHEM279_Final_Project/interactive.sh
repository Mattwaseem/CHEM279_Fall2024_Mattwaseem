#!/bin/sh

docker run --rm -it -v $(pwd):/repo --publish 3000:3000 msse/chem279final
