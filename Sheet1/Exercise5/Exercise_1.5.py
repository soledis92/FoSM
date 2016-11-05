
# coding: utf-8

# In[1]:

import numpy as np


# In[2]:

f = open('numbers.dat', 'rb') # opens numbers file 
data = np.fromfile(f, dtype = float) # read binary data out of file 
data128 = np.fromfile(f, dtype = np.float128) # now with long double data type


# In[3]:

len(data) #make sure it's a million numbers


# In[4]:

### Part a, sum from beginning to end ###
sumup = 0.
for i in range(len(data)):
    sumup =+ data[i]
print('Sum (a) = ', sumup)


# In[5]:

### Part b, reversed direction of summing ###
sumdown = 0.
for i in reversed(range(len(data))):
    sumdown =+ data[i]
print('Sum (b) = ', sumdown)


# In[6]:

### Part c, sort by magnitude (does that mean most negative or the smallest absolute value?) and sum up ###

# Sort data
datasort = np.sort(data) # sort from most negative to most positive
datasortabs = np.sort(np.abs(data)) # sort from smallest absolute to highest absolute

# Check data types 
print('Datasort type: ',datasort.dtype, '\nDatasort absolute data type: ', datasortabs.dtype)

# Sum
sumsortabs = 0.
sumsort = 0.
for i in range(len(datasort)):
    sumsort =+ datasort[i]
    sumsortabs =+ datasortabs[i]
print('Sum (c) =', sumsort, '\nSum (c) absolute =', sumsortabs)


# In[7]:

print(np.sum(np.isnan(data)))


# There are 526 nan values. In both sorted arrays the nan values are at the end of the array. Operations with these values will again give nan results. Therefore both sorted sums are nan. The nan values could be below the double precision. 

# In[8]:

print(np.min(data[np.isnan(data) ==False]), np.max(data[np.isnan(data) == False]))


# The smallest and biggest numbers which are not nan are near to the smallest representable numbers in double precision. A higher precision could avoid the nan values.

# In[9]:

### Part d, as c only with long double data ###

# Sort data
datasort128 = np.sort(data128) # sort from most negative to most positive
datasortabs128 = np.sort(np.abs(data128)) # sort from smallest absolute to highest absolute

# Check data types 
print('Datasort type: ',datasort128.dtype, '\nDatasort absolute data type: ', datasortabs128.dtype)

# Sum
sumsortabs128 = 0.
sumsort128 = 0.
for i in range(len(datasort128)):
    sumsort128 =+ datasort128[i]
    sumsortabs128 =+ datasortabs128[i]
print('Sum (d) =', sumsort128, '\nSum (d) absolute =', sumsortabs128)


# In[10]:

print('Sum beginning to end: ', sumup, '\nSum end to beginning: ',       sumdown, '\nSorted sum most neg to most pos: ',sumsort, '\nSorted sum (absolute values): ', sumsortabs,      '\nSorted sum long double: ', sumsort128, '\nSorted sum long double (abs values): ', sumsortabs128)

