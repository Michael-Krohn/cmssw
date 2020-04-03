def calc_probs(failprobs):
   outprobs = [None] * (len(failprobs)+1)
   if len(failprobs)==1:
      outprobs[0] = 1-failprobs[0]
      outprobs[1] = failprobs[0]
      return outprobs
   prob = failprobs.pop(0)
   fprobs = calc_probs(failprobs)
   outprobs[0] = (1-prob)*fprobs[0]
   outprobs[len(fprobs)]=prob*fprobs[len(fprobs)-1]
   for i in range(1,len(fprobs)):
      outprobs[i] = fprobs[i]*(1-prob)+fprobs[i-1]*prob
   return outprobs
         
      
