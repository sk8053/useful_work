# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 00:42:12 2021

@author: seongjoon kang

"""

import numpy as np
import pandas as pd 
import matplotlib.pylab as plt
from sklearn.metrics import mean_squared_error
import time
import random
from sklearn import gaussian_process
from sklearn.gaussian_process.kernels import RBF, Matern, WhiteKernel
from sklearn.gaussian_process.kernels import ConstantKernel, ExpSineSquared,DotProduct
np.random.seed(12)



data = pd.read_csv('jj.csv')
data = np.array(data)
all_x = data[:,0]
y_all = data[:,1]

test_ind1 = (np.arange(1970, 1975) - 1960)*4
test_ind1 = np.arange(min(test_ind1), max(test_ind1))
test_ind2 = (np.arange(1980, 1991) - 1960) *4
test_ind2 = np.arange(min(test_ind2), max(test_ind2))
# indices of test data
test_ind =  np.append(test_ind1, test_ind2)
# indices of train data
train_ind = list (set(all_x) - set(test_ind))
train_ind = np.array(train_ind, dtype= int)
# get train data from the whole data
X_train = all_x[train_ind-1]
beta_inv = .01
X_train = np.array(X_train, dtype = int)
y_train = data[X_train-1][:,1] +np.random.randn(len(X_train))*np.sqrt(beta_inv)#generative_func(X_train)+np.random.randn(len(X_train))*np.sqrt(beta_inv)

plt.figure()
plt.plot(all_x, y_all, label='true function')
plt.plot(X_train, y_train, '.r', label='training sample (noisy)')
plt.title("Totla Data and Train Data")
plt.xlabel('x')
plt.xticks(np.arange(len(all_x), step = 20), np.arange(1960, 1985, step = 5))
plt.grid()
plt.ylabel('y')
plt.legend()
plt.savefig('data.png')

kernel = 4*RBF(length_scale=20) + 11*ExpSineSquared(periodicity= 10,length_scale=1.1)*RBF(length_scale=13) + 3*WhiteKernel()
#kernel = 0.0015*DotProduct() + ExpSineSquared(periodicity= 10,length_scale=1.1) \
#+ 3*RBF(length_scale = 20) + WhiteKernel()
plt.figure()
plt.imshow(kernel(np.array([all_x]).T))
plt.colorbar()
plt.title('kernel pre-fitting')
plt.savefig('kernel_pre_fitting.png')

gp = gaussian_process.GaussianProcessRegressor(kernel=kernel,normalize_y=True, 
                                               alpha=0,  n_restarts_optimizer=200)
gp.fit(X_train.reshape(-1,1), y_train.reshape(-1,1))
print (gp.kernel_)
plt.figure()
plt.imshow(gp.kernel_(np.array([all_x]).T))
plt.colorbar()
plt.title('kernel post-fitting')
plt.savefig('kernel_post-fitting')


plt.figure()
all_x2 = np.append(all_x, np.arange(85, 120))
mus, sigmas = gp.predict(all_x2.reshape(-1,1), return_std=True)
#plt.plot(all_x2, mus[:,0]+ np.sqrt(sigmas),'k',lw = 0.5)
#plt.plot(all_x2, mus[:,0]- np.sqrt(sigmas),'k',lw = 0.5)
plt.plot(all_x2, mus[:,0],'k', label = 'predicted mean')
plt.fill_between(all_x2, y1 = mus[:,0] + np.sqrt(sigmas), y2 =mus[:,0]- np.sqrt(sigmas), 
                 color='r', label = 'confident range of predicted values')
plt.plot(all_x, y_all, 'b', label = 'original data')


plt.xticks(np.arange(0,len(all_x2)+20, step = 20), np.arange(1960, 1995, step = 5))
plt.xlabel ('Year')
plt.ylabel ('Earning')
plt.title ('Johnson&Johnson Data')
plt.legend()
plt.grid()
plt.savefig('final_predition.png')