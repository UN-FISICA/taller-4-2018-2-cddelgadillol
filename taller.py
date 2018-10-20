import scipy.optimize as sp
import numpy as np

class Derivada:

    def __init__(self,f,metodo="adelante",dx=0.001):
									self.f=f
									self.metodo=metodo;
									self.dx=dx;
    def calc(self,x):
					dx=self.dx
					f=self.f
					if self.metodo=="adelante":
						solucion=f(x+dx)-f(x)
						solucion=solucion/dx
					else:
						if self.metodo=="central":
							solucion=f(x+dx/2)-f(x-dx/2)
							solucion=solucion/dx
						else:
							if self.metodo=="extrapolada":
								f12=(f(x+dx/2)-f(x-dx/2))/dx
								f14=2*(f(x+dx/4)-f(x-dx/4))/dx
								solucion=(4*f14-f12)/3
							else:
								if self.metodo=="segunda":
									f12=(f(x+dx)+f(x-dx)-2*f(x))
									f14=dx*dx
									solucion=(f12)/f14
								else:
									return ("error no valido ingreso")
					return (solucion)
class Zeros:
	def __init__(self,f,metodo,error=1e-4,max_iter=100):
									self.f=f
									self.metodo=metodo;
									self.error=error;
									self.max_iter=max_iter
	def zero(self,vi):
					sol=0;
					solucion=0;
					erra=100;
					f=self.f
					error=self.error
					mamax=self.max_iter
					if self.metodo=="newton":
						while erra>error:
							uno= Derivada(f)
							der=uno.calc(vi)
							vi=vi-f(vi)/der
							erra=f(vi)
							if erra<0:
								erra=erra*(-1)
						return (vi)
					else:
						if self.metodo=="bisectriz":
							a=vi[0]
							b=vi[1]
							aux=0;
							if (f(a)<0 and f(b) < 0) or (f(a)>0 and f(b) > 0):
								return ("error intervalos no validos")
							while erra>error:
								
								medio=(a+b)/2
								uno=f(medio)
								erra=uno
								if uno>0 and f(a)>0:
									a=medio;
								else:
									if uno>0 and f(b)>0:
										b=medio
									else:
										if uno<0 and f(b)<0:
											b=medio
										else:
											if uno<0 and f(a)<0:
												a=medio
											else:
												return(medio)
								
								if erra<0:
									erra=erra*(-1)
							return (medio)
						else:
							if self.metodo=="interpolacion":
								a=vi[0]
								b=vi[1]
								aux=0;
								if (f(a)<0 and f(b) < 0) or (f(a)>0 and f(b) > 0):
											return ("error intervalos no validos")
								while erra>error:
									medio=(a+b)/2
									uno=f(medio)
									erra=uno
									if uno>0 and f(a)>0:
										a=medio;
									else:
										if uno>0 and f(b)>0:
											b=medio
										else:
											if uno<0 and f(b)<0:
												b=medio
											else:
												if uno<0 and f(a)<0:
													a=medio
												else:
													return(medio)
									if erra<0:
										erra=erra*(-1)
								return (medio)
							else:
								if self.metodo=="newton-sp":
									sol=sp.newton(f, vi)
								else:
									if self.metodo=="fsolve-sp":
										sol=sp.fsolve(f,vi) 
									else:
										if self.metodo=="brentq-sp":
											sol=sp.brentq(f, vi[0],vi[1])
										else:
											return "error metodo no valido"
					return(sol)
  

if __name__ == "__main__":
	sol4=0
	sol5=0
	sol6=0
	uno=Derivada(np.sin,"adelante")
	dos=Derivada(np.sin,"central")
	tre=Derivada(np.sin,"extrapolada")
	cua=Derivada(np.sin,"segunda")
	der1=uno.calc(0)
	der2=dos.calc(0)
	der3=tre.calc(0)
	der4=cua.calc(0)
	print("las derivadas son:",der1,der2,der3,der4)
	uno1=Zeros(np.sin,"newton")
	dos2=Zeros(np.sin,"bisectriz")
	tre3=Zeros(np.sin,"interpolacion")
	cua4=Zeros(np.sin,"newton-sp")
	cin5=Zeros(np.sin,"fsolve-sp")
	sei6=Zeros(np.sin,"brentq-sp")
	sol1=uno1.zero(3)
	sol2=dos2.zero((3,3.3))
	sol3=tre3.zero((3,3.3))
	sol4=cua4.zero(3)
	sol5=cin5.zero(3) 
	sol6=sei6.zero((3,3.3))
	print("las zeros estan en :",sol1,sol2,sol3,sol4,sol5,sol6)
