#!/bin/bash
if [ ! -d "Data" ]
then
  mkdir Data
fi
cd Data
echo "Please enter your IU Username"
read USER
sftp "$USER"@carbonate.uits.iu.edu << Script
  cd /N/slate/asolivai/Uploads/AutoClustR
  get -r tpData
  quit
Script