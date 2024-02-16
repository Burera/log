
        
import random
import numpy as np
import matplotlib.pyplot as plt
np.random.seed(0)
m = 6
n = 12
K = n - m
EbN0dB = range(0,10)
Eb_N0_dB = 0


H_matrix = np.array([  [0,	0,	1,	1,	0,	0,	1,	0,	0,	0,	0,	0],
                        [0,	1,	1,	0,	1,	0,	0,	1,	0,	0,	0,	0],
                        [0,	1,	1,	1,	0,	0,	0,	0,	1,	0,	0,	0],
                        [1,	1,	1,	0,	1,	0,	0,	0,	0,	1,	0,	0],
                        [0,	0,	0,	0,	1,	1,	0,	0,	0,	0,	1,	0],
                        [1,	0,	1,	1,	0,	1,	0,	0,	0,	0,	0,	1]])
        
data = np.zeros(K, dtype=np.int8)
data_en = np.zeros(n, dtype=np.int8)
data_de = np.zeros(n, dtype=np.int8)
tx_sig = np.zeros(n, dtype=np.int8)
rx_sig = np.zeros(n, dtype=np.float64)
Lp = np.zeros(n, dtype=np.float64)
LP = np.zeros([m, n], dtype=np.float64)
LP_new = np.zeros([m, n], dtype=np.float64)
Beta = np.zeros([m, n], dtype=np.float64)
Alpha = np.zeros([m, n], dtype=np.int8)
LQ = np.zeros([m, n], dtype=np.float64)
LPmin = np.zeros([m, n], dtype=np.float64)
LQij = np.zeros([m, n], dtype=np.float64)
LPcc = np.zeros([m, n], dtype=np.float64)
LPJ = np.zeros(n, dtype=np.float64)
CH = np.zeros(m,dtype=np.int8)
BER = np.zeros(len(EbN0dB))
LQi = np.zeros_like(H_matrix, float)
NumErrs = 0
NumBits = 0

sum_beta = 0
product_alpha = 0

numiteration = 0
flag = 0
max_iterations = 4
error_counts = np.zeros(4)
ber = np.zeros((4,10))


for Eb_N0_dB in range (0,12):
    
    # Convert EbN0dB to normal value
    EbN0 = np.power(10, (Eb_N0_dB / 10.0))
    Eb = EbN0 * K / n

    # Calculate the value of sigma
    sigMa = np.sqrt(1 / (2 * Eb))

    # Reset to Zero
    NumBits = 0
    error = 0
    error_counts = np.zeros(4)
    while (error_counts[3] < 200):
        
        # Generate dataIn
        for j in range(K):
            data[j] = random.randint(0, 1)
            #print(data_de)
            data_en[j] = data[j]
        
        
        # Encoding (add parity check bit)
        for i in range(K, n):
            data_en[i] = 0
            for j in range(K):
                data_en[i] ^= data[j] * H_matrix[(i - K), j]
        #print("date-en" , data_en)        
        # Modulation
        for j in range(n):
            if data_en[j] == 0:
                tx_sig[j] = 1
            else:
                tx_sig[j] = -1
         # Add noise the transmitting data
        for j in range(n):
            rx_sig[j] = tx_sig[j] + random.gauss(0, sigMa)

        # Enstimate received LLR
        for j in range(n):
            Lp[j] = (2 * rx_sig[j]) / (sigMa * sigMa)

        # Setting the received LLR as initial LP[i][j]
        for i in range(m):
            for j in range(n):
                if H_matrix[i, j] == 1:
                    LP[i, j] = Lp[j] * H_matrix[i, j]
                    
        #print(Lp)
        #print(LP)
        
        for iteration in range(max_iterations):
            
             # Calculate the value of Beta and Alpha
            for i in range(m):
                for j in range(n):
                    if H_matrix[i, j] == 1:
                        Beta[i, j] = np.fabs(LP[i, j])
                        Alpha[i, j] = np.sign(LP[i, j])
                        
            # Estimate LQ[i][j]
            for i in range(m):
                for j in range(n):
                    if H_matrix[i, j] == 1:
                     
                        product_alpha = 1
                        for k in range(n):
                            if (k != j) and (H_matrix[i, k] == 1):
                               
                                product_alpha *= Alpha[i, k]
                           
                        LQ[i, j] = product_alpha 
           
            # Estimate LPmin          
            for i in range(m):
                for j in range(n):
                    if H_matrix[i, j] == 1:
                        min_value = np.inf
                        for k in range(n):
                            if k != j and H_matrix[i, k] == 1:
                                min_value = min(min_value, Beta[i, k])
                        LPmin[i, j] = min_value     
                        
            #Estimate LQij           
            LQij =   LQ *  LPmin        
            
            #  Estimating LP_new
            for j in range(n):
                for i in range(m):
                    if H_matrix[i, j] == 1:
                        LQi[i, j] = 0
                        for k in range(m):
                            if (k != i) and (H_matrix[k, j] == 1):
                                LQi[i, j] += LQij[k, j]
                            LP[i, j] = LQi[i, j] + Lp[j] #UPdate LP here
            
             
            #Calculate LQall 
            LQJ = np.zeros(LQij.shape[1])

            for j in range(LQij.shape[1]):
                indices = np.where(H_matrix[:, j] != 0)[0]
                if len(indices) > 0:
                    LQJ[j] = np.sum(LQij[indices, j])             
            
            
            LPJ = LQJ + Lp
            
            #Hard Decision
            for i in range(n):
                if LPJ[i] > 0:
                    data_de[i] = 0
                   
                else:
                    data_de[i] = 1                                
            #print("data_de" ,data_de)

            #Error Calculation             
            error = 0 
            for i in range(n):
                 if data_en[i] != data_de[i]:
                     error += 1        
            
            #print(error , "error")
            error_counts[iteration] += error  
            iteration += 1    
            
        NumBits += 12
        
    
      
    ber[:,Eb_N0_dB] = error_counts / NumBits       
    print(ber)     
           
        