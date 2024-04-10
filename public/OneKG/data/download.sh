#!/bin/sh

for url in `cat urls.txt`; do
    nohup wget $url &
done
