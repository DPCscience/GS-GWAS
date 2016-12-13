#in terminal type : jupyter qtconsole --style monokai,
from limix.ensemble.lmm_forest import  Forest as LMF
from limix.ensemble import lmm_forest_utils as utils
import scipy as sp
import pylab as pl
# Activate inline plotting
%matplotlib inline
import numpy as np
from sklearn.metrics import r2_score
from sklearn import cross_validation, datasets, linear_model
from scipy.stats.stats import pearsonr

cd ~/Google Drive/MultikernelV6
cd Desktop
kernel = np.loadtxt('kinship_GG_C1_DM.csv', delimiter=',')[0:]
kernel.shape
y= np.loadtxt('pheno_GG_C1_DM.csv', delimiter=',')[:,0:1]
y.shape
x= np.loadtxt('pheno_GG_C1_DM.csv', delimiter=',')

pred=np.loadtxt('pheno_GG_C1_DM.csv', delimiter=',')[:,1:713]

#random state to ensure the same samples are taken each time, is equivalent to set.seed
kf_total = cross_validation.KFold(len(pred), n_folds=5, shuffle=True, random_state=None)
for train, test in kf_total:
	print train, '\n', test, '\n\n'

for train_index, test_index in kf_total:
    print("TRAIN:", train_index, "TEST:", test_index)
    X_train, X_test = pred[train_index], pred[test_index]
    y_train, y_test = y[train_index], y[test_index]
    kernel_train = kernel[sp.ix_(train_index, train_index)]
    kernel_test = kernel[sp.ix_(test_index, train_index)]
	


#kf_total = cross_validation.KFold(len(pred), n_folds=10, shuffle=True, random_state=4)



for k, (train, test) in enumerate(kf_total):
	lm_forest=LMF(kernel='iid',n_estimators=500)
	lm_forest.fit(X_train,y_train)
	response_lmf=lm_forest.predict(X_test)
	print pearsonr(response_lmf,np.reshape(y_test,len(y_test)))


for k, (train, test) in enumerate(kf_total):
	lm_forest=LMF(kernel=kernel_train ,n_estimators=500)
	lm_forest.fit(X_train,y_train)
	response_lmf=lm_forest.predict(X_test,kernel_test)
	print pearsonr(response_lmf,np.reshape(y_test,len(y_test)))

n=0
while n<5 :	
	for k, (train, test) in enumerate(kf_total):
		lm_forest=LMF(kernel=kernel_train,n_estimators=500)
		lm_forest.fit(X_train,y_train)
		response_lmf=lm_forest.predict(X_test,kernel_test)
		print pearsonr(response_lmf,np.reshape(y_test,len(y_test)))
        n +=1	
	

	
def RF():
	for k, (train, test) in enumerate(kf_total):
		print name, 'Starting'
		time.sleep(2)
		lm_forest=LMF(kernel=kernel_train,n_estimators=700)
		lm_forest.fit(X_train,y_train)
		response_lmf=lm_forest.predict(X_test,kernel_test)
		print name, 'Exiting'
		print pearsonr(response_lmf,np.reshape(y_test,len(y_test)))
		


if __name__ == '__main__':
    jobs = []
    for i in range(5):
        p = mp.Process(target=RF)
        jobs.append(p)
        p.start()
	
	

processes = [mp.Process(target=rand_string, args=(5, output)) for x in range(4)]
	
	
	
from  rpy2.robjects import r
r.load("path to .rdata file")
sklearn.cross_validation.KFold(n=4, n_folds=2, shuffle=False,
                               random_state=None)
