from dynspec.tbavg import *
from numpy import *
import pickle
import os

def savePickle(spec,fileroot):
  # saves dyn spec to file in order: nu, t, rr, ll, rl, lr
  nu = spec['freqs']
  t = spec['times']
  rr, rl, lr, ll = [spec['data'][:,i,:] for i in range(4)]

  output = open(fileroot + '.dynspec.pickle','wb')
  # b is for binary
  # w opens and overwrites existing file
  
  pickle.dump(nu,output)
  pickle.dump(t,output)
  pickle.dump(rr,output)
  pickle.dump(ll,output)
  pickle.dump(rl,output)
  pickle.dump(lr,output)
  output.close()

  
def saveTxt(spec,fileroot):
  # creates directory [fileroot].dynspec containing nu.dat, t.dat,
  # rr.dat, rl.dat, lr.dat, ll.dat
  nu = spec['freqs']
  t = spec['times']
  
  pol_list = ['rr','rl','lr','ll']
  pol_data = {}
  for i in range(4):
    pol_data[pol_list[i]] = spec['data'][:,i,:]

  savedir = fileroot + '.dynspec'
  os.system('mkdir ' + savedir)
  savetxt(savedir + '/times.dat',t)
  savetxt(savedir + '/freq.dat',nu)
  for pol in pol_list:
    save(savedir + '/' + pol + '.npy', pol_data[pol])

'''
fileroot = raw_input('Input root file name (e.g., tauCet_Ka): ')
if dir().count('spec')==0:
  print 'loading dynamic spectrum'
  spec = dyn_spec(fileroot + '.tbavg.ms')
saveTxt(spec,fileroot)
savePickle(spec,fileroot)
'''
