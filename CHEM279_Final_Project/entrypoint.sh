#!/bin/sh

set -e

umask 000

# Uses npm to install necessary dependencies like the node api for C++ compilation.
cd /repo/app/frontend
npm install

# Check if build directory exists, if it doesn't make it and compile code.
cd /repo/app/backend
if ! [ -d build ]
then
    mkdir build
fi
cd build
rm -rf *
cmake ../
cmake --build .

# Puts user into the frontend folder so they can easily start the program.
cd /repo/app/frontend

# Eventually this will end with just npm start to build and start web server without user interaction.
echo "Run 'npm start' and ctrl + click server address that appears"

bash "$@"
