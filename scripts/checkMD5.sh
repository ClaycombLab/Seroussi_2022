#!/bin/bash

if [ $# -eq 0 ] || [ ! -e $1 ]
then
    echo "Please pass file with path to raw data files as first arg!"
    exit 1 
fi

while read FILE; do
    # check MD5 checksum, exit if MD5 does not match
    echo -ne "Checking MD5 checksum for ${FILE}: ";
    MD5_CURRENT=`md5sum "${FILE}" | cut -d' ' -f1`;
    MD5_RECORD=`cat "${FILE}.md5" | cut -d' ' -f1`;
    if [ "${MD5_CURRENT}" != "${MD5_RECORD}" ]; then
	echo -e "\e[31mFAILED!!!\e[39m";
	exit 1;
    fi
    echo -e "\e[32mSUCCESS!\e[39m";
done <$1
