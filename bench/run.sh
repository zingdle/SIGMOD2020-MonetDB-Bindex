#!/bin/bash

function run {
  query=${1}
  db=${2}
  method=`echo ${db} | awk -F '-' '{print $2}'`
  runs=20

  # warm up
  mclient -d ${db} -i -t performance -l mal -f trash < mal/${method}/${query}.mal &> /dev/null

  echo "running ${query} ${db}"
  for i in $(seq 1 ${runs}); do
    echo "${i}/${runs}"
    mclient -d ${db} -i -t performance -l mal -f trash < mal/${method}/${query}.mal &> ${query}-${db}-${i}.log
  done
}

run q1 tpch-naive-10
run q1 tpch-bindex-10

run q6 tpch-naive-10
run q6 tpch-bindex-10
