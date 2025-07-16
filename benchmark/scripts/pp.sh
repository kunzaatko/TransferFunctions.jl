#!/bin/bash

awk '
  function ltrim(x) { sub(/^[ \t]*/, "", x); return x; }
  s && NF > 1 && $NF == "["  { s=s $0;               next}
  s && NF == 1 && $1 == "]," { print s "],";   s=""; next}
  s && NF == 1 && $1 == "["  { print s;        s=$0; next}
  s && NF == 1 && $1 == "{"  { print s; print; s=""; next}
  s && NF == 1 && $1 == "]"  { print s $1;     s=""; next}
  s && NF == 1 && $1 == "}"  { print s;        s=$0; next}
  s                          { s=s ltrim($0);        next}
  $NF == "["                 { s=$0;                 next}
  {print}
'
