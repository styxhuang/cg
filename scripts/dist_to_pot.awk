BEGIN { kBT=2.494339 }
NF>=2 {
  g = $2 + 1e-8
  V = -kBT * log(g)
  print $1, V
}