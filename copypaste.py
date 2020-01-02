import numpy as np
from sympy import fwht
import random
import fast_arith as fa
import itertools
import sys
sys.setrecursionlimit(10**6)

class gf2r:
	irr=[]
	r=None
	@staticmethod
	def set_r(x,poly):
		gf2r.r=x
		gf2r.irr=poly #please give poly as  a list of numbers- for eg x^10 + x^2 + 1 is [10,2,1] in any order

	def reduce(self):
		i=None
		r=gf2r.r
		for i in reversed(range(len(self.val))):
			if(self.val[i]==1):
				break
		if i<r and len(self.val)<r:
			self.val=np.append(self.val,np.zeros(r-len(self.val),dtype=np.int64))
			return
		elif i<r and len(self.val)>=r:
			self.val=self.val[0:r]
			return
		while i>=r:
			for x in gf2r.irr:
				self.val[x+i-r]=(self.val[x+i-r]+1)%2
			while self.val[i]==0:
				i=i-1
		self.val=self.val[0:r]

	def reduce_and_return(x):
		i=None
		r=gf2r.r
		for i in reversed(range(len(x))):
			if(x[i]==1):
				break
		if i<r and len(x)<r:
			x=np.append(x,np.zeros(r-len(x),dtype=np.int64))
			return x
		elif i<r and len(x)>=r:
			x=x[0:r]
			return x
		while i>=r:
			#print('lo')
			for u in gf2r.irr:
				x[u+i-r]=(x[u+i-r]+1)%2
			while x[i]==0:
				i=i-1
		x=x[0:r]
		return x

	def __init__(self,vec):
		self.val=np.array(vec)
		self.reduce()

	def __add__(self,other):
		ans=gf2r(list(np.bitwise_xor(self.val,other.val).astype(np.int64)))
		return ans
	def __mul__(self,other):
		a=np.fft.fft(np.append(self.val,np.zeros(gf2r.r-1)))
		b=np.fft.fft(np.append(other.val,np.zeros(gf2r.r-1)))
		c=np.fft.ifft(a*b)
		#print(c)
		c=np.real(c)
		c=c.round()
		c=c.astype(np.int64)
		c=c%2
		ans=gf2r(list(c))
		return ans
		#return gf2r.reduce_and_return(c)
	def __floordiv__(self,other): # p1/p2
		p1=list(self.val)
		p1.reverse()
		p1=list(itertools.dropwhile(lambda x: x==0,p1))
		if(p1==[]):
			p1=[0]
		p2=list(other.val)
		p2.reverse()
		p2=list(itertools.dropwhile(lambda x: x==0,p2))
		if(p2==[]):
			p2=[0]
		elif(p2==[1]):
			ans=gf2r(list(self.val))
			return ans
		# print(p1)
		# print(p2)sudo systemctl start openvpn@iitbvpn.service
		# print('&&')
		q,r=fa.poly_divmod(p1,p2)
		#print(q)
		q1=[int(ab)%2 for ab in q ]
		q1.reverse()
		#print(r)
		ans=None
		if(len(q1)==0):
			ans=gf2r([0])
		else:
			ans=gf2r(q1)
		return ans

	def __mod__(self,other): # p1%p2
		p1=list(self.val)
		p1.reverse()
		p1=list(itertools.dropwhile(lambda x: x==0,p1))
		if(p1==[]):
			p1=[0]
		p2=list(other.val)
		p2.reverse()
		p2=list(itertools.dropwhile(lambda x: x==0,p2))
		if(p2==[]):
			p2=[0]
		elif(p2==[1]):
			ans=gf2r([0]*gf2r.r)
			return ans
		# print(p1)
		# print(p2)sudo systemctl start openvpn@iitbvpn.service
		# print('&&')
		q,r=fa.poly_divmod(p1,p2)
		#print(q)
		r1=[int(ab)%2 for ab in r ]
		r1.reverse()
		#print(r)
		ans=None
		if(len(r1)==0):
			ans=gf2r([0])
		else:
			ans=gf2r(r1)
		return ans


	def power(self,p):
		if p==0:
			x=[0]*gf2r.r
			x[0]=1
			ans=gf2r(x)
			return ans
		if p%2==0:
			a=self.power(p/2)
			return a*a
		else:
			a=self.power((p-1)/2)
			return (a*a)*self
	def inv(self):
		return self.power(int(pow(2,gf2r.r))-2)

	def __truediv__(self,other):
		return self*(other.inv())


	def __str__(self):
		return str(self.val)


def get_bit(x):
	r='{0:b}'.format(x)	
	f=[ int(c) for c in r]
	f.reverse()
	return f

def W(j,x): #finds W_j(x)
	# print(j)
	# print(x)
	ans=gf2r([1])
	for i in range(pow(2,j)):
		y=gf2r(get_bit(i))
		# print(ans)
		# print(x+y)
		# print('#')
		ans=ans*(x+y)
		# print(ans)
	return ans

class fast_transform:
	def __init__(self,h,l): #h is degree of poly and power of 2
		self.h=h
		self.l=l
		self.W2i=[W(i,gf2r(get_bit(int(pow(2,i))))) for i in range(int(np.log2(h)))] #contains W(i,2^i) for i in 0,1,2,... log2h elements
		self.Wil=[W(i,gf2r(get_bit(l))) for i in range(int(np.log2(h)))] #contains W(i,l) for i in 0,1,2,... log2h elements
		self.Wib=[[W(i,gf2r(get_bit(int(b*pow(2,i))))) for b in range(0,int(h/(int(pow(2,i)))),2)] for i in range(int(np.log2(h)))] #2d array first index is i second is b 
		#to give W(i,omega_(b.2^i)), it has logh +1 rows and h/2^(i+1) columns in row i
		#print(self.Wib[0])
		lg=int(np.log2(h))
		print("hello1")
		self.dp=[[[None]*(int(h/int(pow(2,i)))) for m in range(int(pow(2,i)))] for i in range(lg+1)]
		print("hello2")
		#print(dp[0])
	def ft(self,coeffs): #list of coeffs
		lg=int(np.log2(self.h))
		for i in range(self.h):
			self.dp[lg][i][0]=coeffs[i]
		for i in reversed(range(lg)):
			for m in range(int(pow(2,i))):
				for b in range(0,int(self.h/int(pow(2,i))),2):
					
					self.dp[i][m][b]=self.dp[i+1][m][int(b/2)]+\
					((self.Wib[i][b//2]\
						+self.Wil[i])/\
					(self.W2i[i]))\
					*self.dp[i+1][m+int(pow(2,i))][int(b/2)]
				for b in range(1,int(self.h/int(pow(2,i))),2):
					self.dp[i][m][b]=self.dp[i][m][b-1]+self.dp[i+1][m+int(pow(2,i))][b//2]
		return self.dp[0][0]
	def inverse_ft(self,vals):
		lg=int(np.log2(self.h))
		self.dp[0][0]=vals
		for i in (range(lg)):
			for m in range(int(pow(2,i))):
				for b in range(int(self.h/int(pow(2,i+1)))):
					self.dp[i+1][m+int(pow(2,i))][b]=self.dp[i][m][2*b]+self.dp[i][m][2*b+1]
					self.dp[i+1][m][b]=self.dp[i][m][2*b]+((self.Wib[i][b]+self.Wil[i])/self.W2i[i])*self.dp[i+1][m+int(pow(2,i))][b]
		ans=[]
		for m in range(self.h):
			ans.append(self.dp[lg][m][0])
		return ans


def generateAllBinaryStrings(n, i):
	if i == n-1: 
		return [[1],[0]]  
    # First assign "0" at ith position  
    # and try for all other permutations  
    # for remaining positions  
	t1=generateAllBinaryStrings(n,i + 1)
	t2=[x+[1] for x in t1]
	t3=[x+[0] for x in t1]
    # And then assign "1" at ith position  
    # and try for all other permutations  
    # for remaining positions
	return t2+t3
def find_gen(r,poly,id):
	i=1
	while(i<r):
		cand=gf2r(list(map(int,bin(i)[2:])))
		j=1
		fl=1
		while(j<r-1):
			if(np.array_equal(cand.power(j).val,id.val)):
				fl=0
				break
			j=j+1
		if(fl==1):
			return cand
			break
		i=i+1
def list_to_num(l):
	s=0
	po=1
	for i in range(len(l)):
		s=s+l[i]*po
		po=po*2
	return s
def num_to_list(n):
	return list(map(int,bin(n)[2:]))
		
def walsh(a,b,n):
	c=fwht(a)
	d=fwht(b)
	f=[(c[i]*d[i])%(n-1) for i in range(len(c))]
	ans=fwht(f)
	ans=[ans[i]%(n-1) for i in range(len(ans))]
	return ans
def coding(k,n,secrets):
	# secrets is list of gf2r objects of length k
	r1=fast_transform(k,0)
	coeff=r1.inverse_ft(secrets)
	shares=[]
	blocks=int(n/k-1)
	for i in range(blocks):
		r2=fast_transform(k,(i+1)*k)
		newblock=r2.ft(coeff)
		shares=shares+newblock
	# tmp=shares[32:64]
	# r2=fast_transform(k,2*k)
	# coeff2=r2.inverse_ft(tmp)
	# l=[np.array_equal(coeff[i].val,coeff2[i].val) for i in range(k)]
	# print(l)
	return shares

r=10
n=int(pow(2,r))
gf2r.set_r(r,[10,3,0])
# find generator
# l=[0]*r
# ident=gf2r([1])
# g=find_gen(n,gf2r.irr,ident)
# pog={}
# lg={}
# lg[0]=0
# for i in range(n-1):
# 	pog[i]=list_to_num(g.power(i).val)
# 	lg[list_to_num(g.power(i).val)]=i
# print(pog)
# print(lg)

# L=[1,2,1,4,3,3,3,2]
# R=[2,1,3,1,2,1,3,4]
# print(walsh(L,R,8))
# fin=[0,0,0,0,0,0,0,0]
# for i in range(8):
# 	for j in range(8):
# 		fin[i]=(fin[i]+R[j]*L[j^i])%7
# print(fin)
# print(pog)
# print(lg)
# print(len(lg))
# print(lg[1])
# coding
k=32
data=[[random.randint(0,1) for i in range(r)] for j in range(k)]
sec=[gf2r(x) for x in data]
r1=fast_transform(k,0)
shares=coding(k,n,sec)
print(len(shares))
tes2=shares[0:32]

# h=16
# b=gf2r([1,0,1,1,1,0,0,0,1])
# c=gf2r([0,0,0,0,0,0,1,1,1,1,0])
# d=gf2r([1,1,1,0,1,1])
# e=gf2r([0])
# f=gf2r([1,1,1,1,0,1,1])
# g=gf2r([1,1,1,1,1,1,1,0,0,0,1,1,1,1])
# h=gf2r([0,1,0,1,0,1,0,1,1,0,0,1,1,0])
# fr=r.ft([a,b,c,d,e,f,g,h])
# data=[[random.randint(0,1) for i in range(r)] for j in range(h)]
# coeff=[gf2r(x) for x in data]
# r1=fast_transform(h,0)
# pv=r1.ft(coeff)
# print(pv)
# cb=r1.inverse_ft(pv)
# a=[list(coeff[i].val==cb[i].val) for i in range(h)]
# fin=True
# for x in a:
# 	for y in x:
# 		fin = fin and y
# print(fin)


# for x in fr :
# 	print(x)
# print('^')
# bw=r.inverse_ft(fr)
# for x in bw :
# 	print(x)

# data=generateAllBinaryStrings(r,0)
# data=data[:-1]
# coeff=[gf2r(x) for x in data]
# iv=[x.inv() for x in coeff]
# on=[coeff[i]*iv[i] for i in range(len(iv))]
# x=[0]*gf2r.r
# x[0]=1
# ans=gf2r(x)
# fi=[list(x.val==ans.val) for x in on]
# fin=True 
# for i in range(len(fi)):
# 	for j in range(len(fi[i])):
# 		fin=fin and fi[i][j]
# print(fin)








