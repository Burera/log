import random 
import numpy as np
import pdb
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from math import *

M = 2 #BPSK
Eb_N0_dB = 0
num_of_sym = 0
symbol_cx = complex(0,0)
BER_arr = []
symbol_en_real = []
symbol_en_imag = []
symbol_de= np.zeros(9)
n_arr = []
max_iterations = 5 
errors = []
error_counts = np.zeros(5)
ber = np.zeros((5,10))
Lp = np.zeros(9)
m = 4
n = 9
k = n - m
tx_signal = np.zeros(n)
Rx_signal = np.zeros(n)
Lp = np.zeros(n)
LPJ = np.zeros(n)
H = np.array([  [1,	1,	1,	0,	1,	1,	0,	0,	0],
                [1,	1,	1,	1, 	0,	0,	1,	0,	0],
                [1,	1,	0,	1,	1,	0,	0,	1,	0],   
                [1,	0,	1,	1,	1,	0,	0,	0,	1]])

LP = np.zeros([m , n] , float)
LP_new = np.zeros_like(H , float)
Beta = np.zeros_like(H , float)
Alpha = np.zeros_like(H , float)
LQ = np.zeros_like(H , float)
iteration = 5
sum_beta = 0
product_alpha = 0
num_iteration = 0

def  phi(x):
    if x == 0:
        return 0  # Set the result to 0 when x is 0 to avoid division by zero
    else:
        return np.log((np.exp(x) + 1) / (np.exp(x) - 1))
        #return np.log((np.power(np.e , x) + 1) / (np.power(np.e , x) - 1))
# Create a vectorized version of the phi function
vectorized_phi = np.vectorize(phi)
for Eb_N0_dB in range(0,10):
    print('**********************************************************************************************************')
    print("Eb_N0_db : " , Eb_N0_dB)
    
    #Convert Ebno to normal Value
    #Eb_N0_ratio = 10**(Eb_N0_dB/10)*5/9
    Eb_N0_ratio = 10**(Eb_N0_dB/10)
    Eb = Eb_N0_ratio * k /n
    sigma = np.sqrt(1/(2*Eb))

    num_of_sym = 0
    error = 0
    error_counts = np.zeros(5)
    
    while(error_counts[4] < 200 ):
        
        message_bits = np.random.randint(2, size=5)
        
        # multiply message_bits with the first 5 columns of H
        parity_bits = np.mod(np.dot( H[:, :5] , message_bits  ), 2)
        
        # concatenate message_bits and parity_bits to get the codeword
        symbol_en = np.concatenate((message_bits, parity_bits))
        
        modulated_codeword_real = []
        modulated_codeword_imag = []
        
        for bit in symbol_en:
            if bit == 0:
                symbol_en_real = cos(0)
                symbol_en_imag = sin(0)
            elif bit == 1:
                symbol_en_real = cos(pi)
                symbol_en_imag = sin(pi)
            modulated_codeword_real.append(symbol_en_real)
            modulated_codeword_imag.append(symbol_en_imag)
        
        
        #Modulated Signal    
        tx_signal = np.array(modulated_codeword_real)
        #modulated_codeword = np.reshape(modulated_codeword, (1, -1))
        
        #AWGN For LDPC LOG - domain 
        n_arr = []
        for m in range(9):
            AWGN = np.random.normal(0,sigma)
            n_AWGN = AWGN
            n_arr.append(n_AWGN)
            
        Noise = np.array(n_arr)
        
        
        # Add AWGN noise to the modulated codeword
        Rx_signal = tx_signal + Noise
 
        
        
        #Calculate LLR Received
        for i in range(n):
            Lp[i] = (2*Rx_signal[i]) / (sigma * sigma)
            
        #Inserted Received LLR in H matrix
        # Loop through each row of H
        for i in range(H.shape[0]):
            # Find the indices where H equals 1 for the current row
            indices = np.where(H[i] == 1)[0]
            # Set the probability values for the indices where H equals 1 for the current row
            LP[i, indices] = Lp[indices]
            
        #LP = LP_new    
       
        #LP = LP_new
        for num_iteration in range(max_iterations):
            m, n = H.shape
             #Calculate Beta
            for i in range(H.shape[0]):
                #Find the indices where h equals q 
                indices = np.where(H[i] == 1)[0]
                Beta[i , indices] = np.fabs(LP[i, indices])
                
            #Calculate Aplha
            for i in range(H.shape[0]):
                indices = np.where(H[i] == 1)[0]
                Alpha[i , indices] = np.sign(LP[i, indices]) 
                
            #Phi_Beta = vectorized_phi(Beta)
            
            #Product of ALpha
            Product_Alpha = np.zeros_like(Alpha)
        
            # Loop through each row of delta matrix
            for i in range(Alpha.shape[0]):
            # Find the indices where delta does not equal zero for the current row
                indices = np.where(Alpha[i] != 0)[0]
                # Calculate the product of delta values for the current row, if there are nonzero values
                if len(indices) > 0:
                    product = np.prod(Alpha[i, indices])
                    # Set the value of Q for the index of the current row to the product
                    Product_Alpha[i, indices] = product / Alpha[i, indices] 
            
            #summation of beta
            Summation_phi_Beta = np.zeros_like(H, float)
            #aphall = np.zeros_like(H, float)
            for i in range(m):
                for j in range(n):
                    if H[i, j] == 1:
                        sum_beta = 0
                        for k in range(n):
                            if (k != j) and (H[i, k] == 1):
                                sum_beta += phi(Beta[i, k])
                                
                        Summation_phi_Beta[i, j] = (sum_beta)
                        
            Phi_Beta = vectorized_phi(Summation_phi_Beta)
            
            #Calculation of LQij            
            LQij = np.zeros_like(H, float)
            LQij = Product_Alpha * Phi_Beta
            
            
            
            #Calculte LQi
            #  Estimating LP_new[i][j]
            LQi = np.zeros_like(H, float)

            for j in range(n):
                for i in range(m):
                    if H[i, j] == 1:
                        LQi[i, j] = 0
                        for k in range(m):
                            if (k != i) and (H[k, j] == 1):
                                LQi[i, j] += LQij[k, j]

   
            # #Updationg LPij
            #LP = Lp + LQi
            LP = np.zeros_like(H, float)
            for i in range(H.shape[0]):
                for j in range(H.shape[1]):
                    if H[i, j] != 0:
                        LP[i, j] = Lp[j] + LQi[i, j]
            
            #Calculate LQall 
            LQ_all = np.zeros(LQij.shape[1])

            for j in range(LQij.shape[1]):
                indices = np.where(H[:, j] != 0)[0]
                if len(indices) > 0:
                    LQ_all[j] = np.sum(LQij[indices, j])
            
            
            LQJ = LQ_all + Lp
            
            #Hard Decision
            for i in range(n):
                if LQJ[i] > 0:
                    symbol_de[i] = 0
                   
                else:
                    symbol_de[i] = 1
                    
                    
            #Error Calculation             
            error = 0 
            for i in range(n):
                 if symbol_de[i] != symbol_en[i]:
                     error += 1        
            
            error_counts[num_iteration] += error  
            num_iteration += 1
            
            
        num_of_sym += 9    
            
    ber[:,Eb_N0_dB ] = error_counts / num_of_sym   
       