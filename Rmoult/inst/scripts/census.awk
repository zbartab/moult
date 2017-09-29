#!/usr/bin/awk -f

BEGIN {
  y = -1;
}

{
  if ($2 == 51) { 
    if ($14 == 0) census0[$1]++;
    else if ($14 == 1) census1[$1]++;
    else census2[$1]++;
  }
  if ($4 ~ "PREDA") {
    if ($14 == 0) pred0[$1]++;
    else if ($14 == 1) pred1[$1]++;
    else pred2[$1]++;
  }
  if ($4 ~ "STARV") {
    if ($14 == 0) starv0[$1]++;
    else if ($14 == 1) starv1[$1]++;
    else starv2[$1]++;
  }
}

END {
  for (year in census2) {
    printf "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", 
      year, census0[year], census1[year], census2[year], pred0[year], 
      pred1[year], pred2[year], starv0[year], starv1[year],starv2[year];
  }
}
