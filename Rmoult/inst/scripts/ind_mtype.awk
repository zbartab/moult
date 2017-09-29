#!/usr/bin/awk -f

BEGIN {
  if (y == "") {
    y = 2;
  }
}

{
  if ($1 == y) {
#    if ($7 ~ "START") {
#      sched[$3] = sprintf("%s%s", sched[$3], "B");
#    }
    if ($7 ~ "START") {
      sched[$3] = sprintf("%s%s", sched[$3], "B");
    } else if ($7 ~ "ABAND") {
      sched[$3] = sprintf("%s%s", sched[$3], "a");
    } else if ($7 ~ "ABORT") {
      sched[$3] = sprintf("%s%s", sched[$3], "d");
    }
    if ($5 ~ "START") sched[$3] = sprintf("%s%s", sched[$3], "p");
    if ($6 ~ "START") sched[$3] = sprintf("%s%s", sched[$3], "s");
    if ($2 == 51) surv[$3]++;
  } else if ($1>y && FNR > 1) exit;
}

END {
  for (bird in surv) {
    if(sched[bird] == "") printf "%s\t%s\n", bird, "n";
    else printf "%s\t%s\n", bird, sched[bird];
  }
}
