
# coding: utf-8

# ## Elastic Properties
# ### Code for extracting elastic tensor and calculating mechanical properties from VASP OUTCAR
# 
# Equations can be found at https://www.materialsproject.org/wiki/index.php/Elasticity_calculations

# In[86]:

import numpy as np


# In[87]:

def get_elastic_tensor(filename):
    ''' Reads the elastic tensor from the OUTCAR. 
    Args:
        filename : the name of the vasp OUTCAR
    Returns:
        elastic_tensor : 6x6 tensor of the elastic moduli
    '''
    f = open(filename,"r")
    lines = f.readlines()
    f.close()
    copy = False
    elastic_tensor = []
    for line in lines:
        inp = line.split()
        if inp == []:
            continue 
        if len(inp) < 4 or len(inp) > 7:
            continue
        if len(inp) == 4 and inp[0] == 'TOTAL':
            copy = True
        if copy:
            if len(inp) == 7 and len(inp[0]) == 2:
                elastic_tensor.append(inp[1:])
    return np.asarray(elastic_tensor).astype(np.float)


# ### Elastic tensor $C_{ij}$

# In[88]:

elastic_tensor = get_elastic_tensor('OUTCAR')


# ### Divide by 10 to convert kBar to GPa

# In[89]:

Cij = elastic_tensor/10

# ### Compliance tensor $s_{ij}$ $(GPa^{-1})$
# $s_{ij} = C_{ij}^{-1}$

# In[90]:

Sij = np.linalg.inv(Cij)


# ### Voigt bulk modulus $K_v$ $(GPa)$
# $9K_v = (C_{11}+C_{22}+C_{33}) + 2(C_{12} + C_{23} + C_{31}) $

# In[91]:

Kv = ((Cij[0,0] + Cij[1,1] + Cij[2,2]) + 2 * (Cij[0,1] + Cij[1,2] + Cij[2,0])) / 9


# ### Reuss bulk modulus $K_R$ $(GPa)$
# $1/K_R = (s_{11}+s_{22}+s_{33}) + 2(s_{12} + s_{23} + s_{31})$

# In[92]:

Kr = 1/((Sij[0,0] + Sij[1,1] + Sij[2,2]) + 2 * (Sij[0,1] + Sij[1,2] + Sij[2,0])) 


# ### Voigt shear modulus $G_v$ $(GPa)$
# $15 G_v = (C_{11}+C_{22}+C_{33}) - (C_{12} + C_{23} + C_{31}) + 3(C_{44} + C_{55} + C_{66})$

# In[93]:

Gv = (4 * (Cij[0,0] + Cij[1,1] + Cij[2,2]) - 4 * (Cij[0,1] + Cij[1,2] + Cij[2,0]) + 3 * (Cij[3,3] + Cij[4,4] + Cij[5,5]))/15


# ### Reuss shear modulus $G_v$ $(GPa)$
# $ 15/G_R = 4(s_{11}+s_{22}+s_{33}) - 4(s_{12} + s_{23} + s_{31}) + 3(s_{44} + s_{55} + s_{66})$

# In[94]:

Gr = 15 / (4 * (Sij[0,0] + Sij[1,1] + Sij[2,2]) - 4 * (Sij[0,1] + Sij[1,2] + Sij[2,0]) + 3 * (Sij[3,3] + Sij[4,4] + Sij[5,5]))


# ### Voigt-Reuss-Hill bulk modulus $K_{VRH}$ $(GPa)$
# $K_{VRH} = (K_R + K_v)/2$

# In[95]:

Kvrh = (Kv + Kr)/2
print "Voigt-Reuss-Hill bulk modulus (GPa): ", Kvrh


# ### Voigt-Reuss-Hill shear modulus $G_{VRH}$ $(GPa)$
# $G_{VRH} = (G_R + G_v)/2$

# In[96]:

Gvrh = (Gv + Gr)/2
print "Voigt-Reuss-Hill shear modulus (GPa): ", Gvrh


# ### Isotropic Poisson ratio $\mu$
# $\mu = (3K_{VRH} - 2G_{VRH})/(6K_{VRH} + 2G_{VRH})$

# In[83]:

mu = (3 * Kvrh - 2 * Gvrh) / (6 * Kvrh + 2 * Gvrh )
print "Isotropic Poisson ratio: ", mu


# In[85]:

Cij


# In[ ]:



