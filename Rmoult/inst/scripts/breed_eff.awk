#!/usr/bin/awk -f

BEGIN {
	if(y=="") y = 2;
}

{
  if ($1==y-1 && $2==51 && $4 ~ /^ALIVE/ && $14==2) cp[$3]++;
  if ($1==y && $2==51 && $4 ~ /^ALIVE/ && $14==2 && cp[$3]>0) c[$3]++;
  if ($1==y && ($7 ~ /^START/ || $7 ~ /^KEEP/) && cp[$3]>0) eff[$3]++;
  if ($1==y && $7 ~ /^START/ && cp[$3]>0) bs[$3]++;
  if ($1==y && $7 ~ /^ABAND/ && cp[$3]>0) ba[$3]++;
}

END {
  c_b = eff_b = bs_b = ba_b = 0;
  for(bird in c) {
    c_b += c[bird];
    eff_b += eff[bird];
    bs_b += bs[bird];
    ba_b += ba[bird];
  }
  if(bs_b > 0) ba_b /= bs_b;
  if(c_b > 0) {
    eff_b /= c_b;
    bs_b /= c_b;
  }
  print c_b, eff_b, bs_b, ba_b;
}

