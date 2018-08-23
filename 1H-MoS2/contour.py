from pylab import *
from math import *
nx = 300
ny = 300
index='BERRYCURV'
state='233'
#if len(sys.argv) == 2:
# state=sys.argv[1]
#elif len(sys.argv) == 3:
# index=str(sys.argv[1])
# state=str(sys.argv[2])

name=index+state
#filein=index+'.'+state+'.dat'
#fileou=index+'.'+state+'.eps'
filein=index+'.dat'
fileou=index+'.eps'
#if sys.argv[1] == '-cd':
# filein='CIRC_DICHROISM.dat'
# fileou='CIRC_DICHROISM.eps'
#print '# INPUT  LIST :', filein
data = np.loadtxt(filein)
colorlevel = 51
x,y,z = data[:,0], data[:,1], data[:,3]

# If you want log10 scale plot 
#for i in range(len(z)):
# z[i-1]=(log10(abs(data[i-1,3])))

#data_k = np.loadtxt('./SKP.dat')
#x1,x2 = data_k[:,0], data_k[:,1]
#data_kgrd = np.loadtxt('./1.dat')
#g1,g2 = data_kgrd[:,0], data_kgrd[:,1]

xi = np.linspace(x.min(),x.max(),nx)
yi = np.linspace(y.min(),y.max(),ny)
zi = griddata(x,y,z,xi,yi)
maxz=z.max()
minz=z.min()
deltaz=(maxz-minz)/colorlevel
level=arange(minz,maxz+deltaz/5,deltaz)
#level=arange(-16.0,16.0,0.1)
#level_1=arange(-0.004,-0.0003,0.0004)
#level_2=arange(0.0004,0.0041,0.0004)
#level_3=arange(-0.0004,0.00041,0.0008)
#CS= plt.contour(xi,yi,zi,level,linewidths=0.5,colors='k',extend='both')
#CS= plt.contourf(xi,yi,zi,level_3,alpha=0.75)
#CS= plt.contour(xi,yi,zi,level_1,linewidths=0.25,colors='k',extend='both')
#CS= plt.contour(xi,yi,zi,level_2,linewidths=0.25,colors='k',extend='both')
#CS= plt.contourf(xi,yi,zi,level,alpha=0.75,cmap=cm.RdYlBu)
CS= plt.contourf(xi,yi,zi,level,alpha=0.75,cmap=cm.coolwarm)
#CS= plt.contourf(xi,yi,zi,level_3,alpha=0.75,colors='White')
#plt.clabel(CS,colors='black', inline=1,fontsize=5)

bounds = np.linspace(z.min(),-1*z.min(),10)
#bounds = [-70,-50,-30,-10,10,30,50,70]
plt.colorbar(ticks=bounds,format='%5.2f')

#scatter(x1,x2,s=10,alpha=0.75,zorder=2)
scatter(x,y,s=.2,alpha=0.75,zorder=2)
xin=x.min()#+2.500969/2
xen=x.max()#-2.500969/2
yin=y.min()#+2.887871/1
yen=y.max()#-2.887871/1
plt.axis([xin,xen,yin,yen])
plt.setp(gca(),yticks=(arange(0)),xticks=(arange(0)),aspect='equal')
setp(gca(), yticklabels=[], yticks=(), xticks=())
#clabel(CS,inline=1, fontsize=10)

#savefig('./contour.png',dpi=96)
savefig(fileou,bbox_inches='tight',transparent=True,pad_inches=0)
print '# OUTPUT LIST :', fileou
