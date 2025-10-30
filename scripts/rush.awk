BEGIN {
  # Construct temporary files
  cmd = "mktemp"
  cmd | getline
  close(cmd)
  seq = $1
  cmd | getline
  close(cmd)
  query = $1
  cmd | getline
  sbjct = $1
  close(cmd)
  # Construct simulation command
  template = "ms 2 1 -t 1000 -r %f 100000 | ms2dna > %s; "
  template = template "getSeq S1 %s > %s; "
  template = template "getSeq S2 %s > %s; "
  template = template "rush -q %s %s | grep -v '^s'"
  minRho = 0
  maxRho = 4096
  it = 100
  print "#rho\tRejection (alpha = 0.05)"
  # Iterate over values of \rho
  for (r = minRho; r <= maxRho; r *= 2) {
    cmd = sprintf(template, r, seq,
		  seq, query,
		  seq, sbjct,
		  query, sbjct)
    c = 0
    for (j = 0; j < it; j++) {
      cmd | getline
      if ($6 <= 0.05)
	c++
      close(cmd)
    }
    print r "\t" mean "\t" c/it
    if (r == 0)
      r = 0.5
  }
  cmd = "rm " seq " " query " " sbjct
  system(cmd)
}
