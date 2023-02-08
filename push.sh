#!/bin/bash
CurrentDate=`date +"%d/%m"`

git add .
git commit -m $CurrentDate
git push

