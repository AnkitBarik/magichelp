#!/bin/bash

function getCheck(){
   ls -lth check*|head -n 1|awk '{print $9}'
}

function getTag(){
   tagNum=`grep tag input.nml |awk -F'run' '{print $2}'|tr -dc '[0-9]+'`
   tagNumNew=$(($tagNum + 1))
   echo $tagNumNew
}
