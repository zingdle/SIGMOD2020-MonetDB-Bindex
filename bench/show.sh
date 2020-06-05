#!/bin/bash

function show {
  query=${1}
  db=${2}
  method=`echo ${db} | awk -F '-' '{print $2}'`
  time=`cat ${query}-${db}-*.log | grep -oP "(?<=run:)\w+\.\w+" | awk '{ total += $1 } END { print total/NR }'`

  echo ${query} ${method} ${time}
}

show q1 tpch-naive-10
show q1 tpch-bindex-10

show q6 tpch-naive-10
show q6 tpch-bindex-10
