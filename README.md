# Cálculo para la potencias de frenados en agua.
Se calculó la potencia de frenado de protones con energías entre 1-250 MeV, se ajustaron diferentes modelos a los datos obtenidos de PSTAR. De igual forma con una simple rutina se hizo un cálculo aproximado resolviendo la ecuación de Bethe-Bloch.
## Es necesario tener las siguientes bibliotecas de python instaladas 
'''
numpy
matplotlib.pyplot
scipy.optimize
scipy.stats
'''
### Cálculo de la ecuación de Bethe-Bloch sin correcciones
'''
#definimos las constantes a utilizar
k = 0.307075 # MeV g^-1 cm^2
z = 1 # carga del protón en e
ZA = 0.55509 #Z/A 
c = 2.9979e8  #velocidad de la luz
M =  939.595 # masa del protón en MeV/c 
me = 0.511 # masa del electrón
rho = 1 # g cm^-3 densidad del material en este caso agua
I = 7.5e-5 #MeV
E0 = 100
#definimos algunas funciones
def Tmax(beta,gama):
    return (2*me*pow(beta*gama,2))*pow(10,6)/(1+2*gama*(me/M)+pow((me/M),2))
def beta2(E):
    return np.sqrt(pow(M,2)-pow(E,2))/np.sqrt(pow(M,2)+pow(E,2))
#Imtegración númerica de la ec. de  Bethe-Bloch en el rango de energía seleccionado
E0 = 0
Elost = 0
E=E0
dx = 0.01  #paso inicial
x = 0.0    #valor inicial de x
m = 0    #contador
n = 2500
xx = np.zeros(n)  #arrays para guardar los resultados
yy = np.zeros(n) 
print('Energy(MeV)\t\tSP')   #data table header
for E in np.arange(1,250,2):
    x+= dx  #update x-position
    beta = beta2(E)
    gama = 1/(1-pow(beta,2))
    tmax = Tmax(beta,gama)
    #Resolvemos la ec. para los valores de beta y gamma
    dEdx = k*pow(z,2)*ZA*rho*(1/pow(beta,2))*(0.5*np.log(2*0.511*pow(beta*gama,2)*tmax/(pow(I,2)))-pow(beta,2))
    E -= dEdx*dx
    Elost+=dEdx*dx
    xx[m]=x
    yy[m]=dEdx
    m=m+1
    print('%f \t %f'% (E,dEdx))
'''
### Para realizar el ajuste
Cargamos los datos a utilizar:
'''
#Load the data set from PSTAR-NIST
PSTAR = pd.read_csv('PSTAR.txt',sep='\s+',header=None)
PSTAR = pd.DataFrame(PSTAR)
a1 = PSTAR[0] #Energy T -> MeV
a2 = PSTAR[1]/1000.# Stopping Power Elec MeV/cm²g -> convert to mg
a3 = PSTAR[2]  # Nuclear sp
a4 = PSTAR[3]  # Total SP 
a5 = PSTAR[4]  #CSDA Range
a6 = PSTAR[5]  #Projected Range
'''
Despues definimos los modelos y la estimación inicial de parametros.
'''
#Bragg-Klemann Rule
def BG(x, p, a):
    return (x**(1-p))/p*a
#Guessed parameters:
def get_rsq(f, y, popt):
    ss_res = np.dot((a2 - BG(a1, *popt)),(a2 - BG(a1, *popt)))
    ymean = np.mean(a2)
    ss_tot = np.dot((a2-ymean),(a2-ymean))
    return 1-ss_res/ss_tot
p0 = [1, 0.5]    
popt,pcov = op.curve_fit(BG, a1, a2, p0=p0)     
print("Estimated Parameters",popt)  
print("Mean R :",  get_rsq(BG, a2, popt))
#compute one standard deviation errors on the parameters use 
perr = np.sqrt(np.diag(pcov))
print("Std dev errors", perr)
x_BG = np.linspace(0,a1.max(),50)
a,p = popt
y_BG = BG(x_BG,a,p)
'''
Despues gráficamos el resultado
'''
#Plot the results
plt.plot(x_BG,y_BG,'k-', label='Bragg-Kleeman')
plt.plot(a1,a2, '--', label = 'PSTAR Data')
plt.yscale('log')
plt.xlabel('Proton energy E (MeV)')
plt.ylabel('Stopping power SP(E)'' $(MeV cm^{2}/g)$')
#plt.xscale('log')
plt.grid()
plt.legend(loc='best')
#plt.savefig('braggkl.png')
plt.show()
'''
The code in python consist of a simple routine to calculate the stopping power for protons in water. Then it loads the data 
from the data folder, plot them and fit the data to two models for the calculation of stopping powers, the Bethe-Bloch formula and
the Bragg-Kleeman rule.
