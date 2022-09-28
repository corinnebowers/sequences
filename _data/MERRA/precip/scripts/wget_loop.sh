#!/bin/bash

for url in $(cat url_subset.txt) 
do ./script.exp $url
done