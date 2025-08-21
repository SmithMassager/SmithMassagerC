#!/bin/bash

for i in 101 1009 10007 50021 100003 1000003 10000019
do
  time ./src/a.out $i &> $i.out
done
