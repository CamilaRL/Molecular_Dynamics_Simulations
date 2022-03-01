import numpy as np
import matplotlib.pyplot as plt

endereco_E = "C:\\Users\\camil\\Desktop\\UFRGS\\IC-Fluidos\\tButanol\\PROGRAMAS\\Simulacoes\\energias.txt"
U, K, etot= np.loadtxt(endereco_E, skiprows = 1, unpack = True)

endereco_G = "C:\\Users\\camil\\Desktop\\UFRGS\\IC-Fluidos\\tButanol\\PROGRAMAS\\Simulacoes\\g(r).txt"
g, r = np.loadtxt(endereco_G, skiprows = 1, unpack = True)

t = np.arange(0, 1000)

# Gráfico da Energia pelo tempo
plt.plot(t, U, 'r.', label='Energia Potencial (U)')
plt.plot(t, K, 'g.', label='Energia Cinética (K)')
plt.plot(t, etot, 'k.', label='Energia Total (E)')
plt.xlabel('Tempo')
plt.ylabel('Energia')
plt.ylim(bottom=-5, top=5)
plt.legend(loc='right')
plt.title('Energia do Sistema', fontweight='bold')
plt.show()

plt.plot(r, g, 'k.', label='nhis: 15')
plt.xlabel('r')
plt.ylabel('g(r)')
plt.legend()
plt.show()
