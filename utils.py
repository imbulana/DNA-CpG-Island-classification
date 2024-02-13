def cgi_partition(file_path, cpi): # use generator?
  seq=''
  P0=[] # non CpG
  P1=[] # CpG
  with open(file_path, 'r') as f:
    for line in f:
      if not line.startswith('>'):
        seq+=line.strip()

  prev_e=0
  for s,e in cpi:
    if s-prev_e>150: # ignore if <150 bp
      P0.append(seq[prev_e:s])
    P1.append(seq[s:e+1])
    prev_e=e+1
  
  if len(seq)-prev_e>150:
    P0.append(seq[e+1:])
  
  # partition P0 and P1 into sequences of length 150
  X0=[]
  X1=[]
  for s in P0:
    partition(s, X0, 150)
  for s in P1:
    partition(s, X1, 150)
  return X0, X1

def partition(s, P, length):
  j=0
  n_partitions=len(s)//length
  for _ in range(n_partitions):
    P.append(s[j:j+length])
    j+=length

# cgi_partition('AL031718.fas', [(4620, 8266), (16033, 18297)])
