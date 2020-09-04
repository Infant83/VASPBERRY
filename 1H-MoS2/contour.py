from pylab import *
from math import *
from scipy.interpolate import griddata
nx = 300
ny = 300

filein='BERRYCURV.dat'
fileou='BERRYCURV.eps'
data = np.loadtxt(filein)
colorlevel = 51
x,y,z = data[:,0], data[:,1], data[:,3]

xi = np.linspace(x.min(),x.max(),nx)
yi = np.linspace(y.min(),y.max(),ny)
zi = griddata((x,y),z,(xi[None,:],yi[:,None]))
maxz=z.max()
minz=z.min()
deltaz=(maxz-minz)/colorlevel
level=arange(minz,maxz+deltaz/5,deltaz)
#CS= plt.contourf(xi,yi,zi,level,alpha=0.75,cmap=cm.coolwarm)
CS= plt.contourf(xi,yi,zi,level,cmap=cm.coolwarm)

bounds = np.linspace(z.min(),-1*z.min(),10)
plt.colorbar(ticks=bounds,format='%5.2f')

#scatter(x,y,s=.2,alpha=0.75,zorder=2)
scatter(x,y,s=.2,zorder=2)
xin=x.min()#+2.500969/2
xen=x.max()#-2.500969/2
yin=y.min()#+2.887871/1
yen=y.max()#-2.887871/1
plt.axis([xin,xen,yin,yen])
plt.setp(gca(),yticks=(arange(0)),xticks=(arange(0)),aspect='equal')
setp(gca(), yticklabels=[], yticks=(), xticks=())
savefig(fileou,bbox_inches='tight',transparent=True,pad_inches=0)
