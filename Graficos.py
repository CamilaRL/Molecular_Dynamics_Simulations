import numpy as np
import matplotlib.pyplot as plt

endereco = "C:\\Users\\camil\\Desktop\\UFRGS\\IC-Fluidos\\tButanol\\PROGRAMAS\\Simulacoes\\saida.txt"
U, K, etot= np.loadtxt(endereco, skiprows = 1, unpack = True)

t = np.arange(0, 1000)

'''
fig = plt.figure(figsize=(20, 5))

plt.subplot(131)
plt.plot(t, U, 'r.')
plt.xlabel('Tempo')
plt.ylabel('Energia')
plt.title('Energia Potencial')
plt.subplot(132)
plt.plot(t, K, 'g.')
plt.xlabel('Tempo')
plt.ylabel('Energia')
plt.ylim
plt.title('Energia Cinética')
plt.subplot(133)
plt.plot(t, etot, 'k.')
plt.xlabel('Tempo')
plt.ylabel('Energia')
plt.title('Energia Total')

plt.subplots_adjust(wspace=0.35)
plt.suptitle('Energia do sistema', fontweight='bold')
plt.show()
'''

plt.plot(t, U, 'r.', label='Energia Potencial (U)')
plt.plot(t, K, 'g.', label='Energia Cinética (K)')
plt.plot(t, etot, 'k.', label='Energia Total (E)')
plt.xlabel('Tempo')
plt.ylabel('Energia')
plt.ylim(bottom=-5, top=5)
plt.legend(loc='right')
plt.title('Energia do Sistema', fontweight='bold')
plt.show()