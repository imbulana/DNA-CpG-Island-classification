import numpy as np

class Markov:
  def __init__(self):
    '''
    - T0 is the transition matrix for non CpG island model
    - T1 is the transition matrix for CpG island model
    '''
    self.T0=None
    self.T1=None

  def train(self, X0, X1, transitions):
    '''
    Fill transition matrices T0, T1
    - X is a list of DNA sequences
    - transitions is square matrix containing all possible transitions
    '''
    k=len(transitions)
    self.T0=[[0 for _ in range(k)] for _ in range(k)]
    self.T1=[[0 for _ in range(k)] for _ in range(k)]

    def avg(X, t):
      # m=len(t)
      def seq_avg(seq):
        n=len(seq)
        count=seq.count(t)
        # return pseudocount, (count+1)/((n-m+1) transitions in seq + 12 possible transitions)
        return (count+1)/(n+11)
      
      sum_avg=0
      for seq in X:
        sum_avg+=seq_avg(seq)
      return sum_avg/len(X)

    def fill(X, T):
      for i in range(k):
        for j in range(k):
          if i!=j:
            T[i][j]=avg(X, transitions[i][j])
        T[i][i]=1-sum(T[i])
    
    fill(X0, self.T0)
    fill(X1, self.T1)

  def predict(self, seq):
    '''
    - seq is a DNA sequence
    '''
    M={'A':0, 'C':1, 'G':2, 'T':3}
    def l(T):
      log_likelihood=0 # include probability for start state log(pi[seq[0]])?
      for i in range(len(seq)-1):
        log_likelihood+=np.log(T[M[seq[i]]][M[seq[i+1]]])
      return log_likelihood
    
    l0=l(self.T0) # log(P(seq | non CpG))
    l1=l(self.T1) # log(P(seq | CpG))

    return self.binary_classifier(l0, l1)

  def predict_all(self, X):
    '''
    - X is a list of DNA sequences
    '''
    P=np.array([0]*(len(X)))
    for i in range(len(X)):
      P[i]=self.predict(X[i]) 
    return P

  def binary_classifier(self, l0, l1):
    '''
    Classify a sequence as CpG=1 or non CpG=0 island with Naive Bayes
    - l0 is the log likelihood of getting sequence via non CpG model
    - l1 is the log likelihood of getting sequence via CpG model
    '''
    # priors
    lp0=np.log(.75) # log(P(non CpG))
    lp1=np.log(.25) # log(P(CpG))

    x=l1-l0+lp1-lp0 # log(P(seq | CpG)/P(seq | non CpG)) + log(P(CpG)/P(non CpG))
    return 1 if .5+.5*np.tanh(x/2)>.5 else 0 # return 1 if sigmoid(x)>.5 else 0