#!/usr/bin/awk -f

BEGIN {
        c = p = s = 0;
	cy = py = sy = 0;
	if(y=="") y = 2;
}

{
  if ($1 == y-1 && $2 == 51 && $4 ~ /^ALIVE/) {
    cy++;
    if($14 == 2) cye[$3]++;
  }
  if ($1==y && $2==51 && $4 ~ /^ALIVE/) c++;
  if ($1<y+1) {
    if($4 ~ /^PREDA/) p++;
    if($4 ~ /^STARV/) s++;
  }
  if ($1==y) {
    if($4 ~ /^PREDA/) {
      py++;
      if($14 == 2) pye[$3]++;
    }
    if($4 ~ /^STARV/) {
      sy++;
      if($14 == 2) sye[$3]++;
    }
  }
##  if($1 > y && NR > 1) nextfile;
}

END {
  cyes = pyes = syes = 0;
  for(bird in cye) {
    cyes += cye[bird];
    pyes += pye[bird];
    syes += sye[bird];
  }
  print c, p, s, cy, py, sy, cyes, pyes, syes;
}

