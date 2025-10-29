import matplotlib.pyplot as plt
import random
from PIL import Image
import numpy as np
from Isotropic.Metropolis import Metropolis
from Isotropic.Gibbs import Gibbs
from Anisotropic.AMetropolis import AMetropolis
from Anisotropic.AGibbs import AGibbs


#convert an image to a binary matrix debending otn the specified threshold
#RETURN the sample containing the image,  height and width of it without considering the additional borders (for coerence with the rest of the code)
def img_to_matrix(img_path,threshold = 128):
    # Read the image
    img = Image.open(img_path)

    # Convert the image to a NumPy array
    img_array = np.array(img)

    # Threshold the pixel value
    height, width = img_array.shape[:2]
    height = height + 2 #add border
    width = width + 2 #add border
    sample = np.zeros((height)*(width))
    for i in range(1,height-1):
        for j in range(1,width-1):

            if img_array[i-1][j-1][0] > threshold: #here no border
                sample[i*width+j] = 1
            else:
                sample[i*width+j] = -1

    return sample,width-2,height-2


#add a noise to the image with some flipping probability
#Return the BW image with the added noise
def add_noise(matrix, flip_prob):
     # Define a function to flip values randomly with probability flip_prob
        flip_func = np.vectorize(lambda x: -x if np.random.rand() < flip_prob else x)

        # Apply the flip function to each element in the matrix
        noisy_matrix = flip_func(matrix)

        return noisy_matrix

def print_sample(Sample, width, height):
    for i in range(1,height-1):
        for j in range(1,width-1):
            if Sample[j + i * width] == -1:
                print('o', end=' ')
            elif Sample[j + i * width] == 1:
                print('x', end=' ')
            else:
                print(' ', end=' ')
        print()  # Print a newline after each row

# Example usage:




Operation ='Create a sample of the Ising widthRF model'
'''Alpha = float(input('Choose the label cost parameter Alpha\nAlpha=0 : equally probable labels -1 and 1; \nAlpha>0: preference for -1,  Alpha<0: preference for +1 (try ~ 0.05) ;\nAlpha= '))
Beta = float(input('Choose the interaction parameter Beta\nBeta<0 labels of the same type cluster (try Beta=-0.9)\nBeta= '))
width = int(input('Insert the number of nodes in the horizontal direcition \n A good value is 100\n width= '))
height = int(input('Insert the number of nodes in the Vertical direcition \n A good value is 100\n height= '))'''
Alpha = 0 #site field
Beta = -0.8 #interaction
width = 120
height = 100
Beta_h = -0.1
Beta_v = -0.1
Beta_d1 = -0.1
Beta_d2 = -0.8
#Uncomment to have an example with image denoising
img_path = 'Isotr' \
'opic/txt2.jpeg' #path to the test image
Sample,width,height = img_to_matrix(img_path)
Sample = add_noise(Sample, 0.2)

ITERA = 50 #number of total iterations is this*number of sites
T=0.1


L = (width+2)*(height+2) #add border
#Sample = [0] * L
generate_other = 'y'
while generate_other.lower() == 'y':
    '''for i in range(height+2):
        for j in range(width+2):
            u = random.randint(0, 1)
            Sample[j + i * (width+2)] = 2 * u - 1  # 2*u-1 is a trick to convert 0 to -1 and 1 to 1
'''
    Sample2 = Sample.copy() #use the same initial sample also for Gibbs


    Z = Gibbs(Sample,height,width,0,Alpha,Beta,T) #just to obtain Z and show the initial graph
    plt.imshow(Z, cmap='gray')
    plt.show()
    X=Metropolis(Sample,height,width,ITERA,Alpha,Beta,T,lam=0.1, T0=0.5,Tf=0.01,anneal=False)
    #X=AMetropolis(Sample,height,width,ITERA,Alpha,Beta_h,Beta_v,Beta_d1,Beta_d2,T,lam=0.3, T0=0.5,Tf=0.01,anneal=False)
    plt.imshow(X, cmap='gray')
    plt.title(f"Montecarlo beta={Beta},Iterations={ITERA}")
    #plt.title(f"Montecarlo Beta_h={Beta_h},Beta_v={Beta_v},Beta_d1={Beta_d1},Beta_d2={Beta_d2},Iterations={ITERA}")
    plt.savefig(f"Montecarlo_{ITERA}")
    plt.show()
    Y = Gibbs(Sample,height,width,ITERA,Alpha,Beta,T,lam=0.1, T0=0.5,Tf=0.01,anneal=False)
    #Y = AGibbs(Sample,height,width,ITERA,Alpha,Beta_h,Beta_v,Beta_d1,Beta_d2,T,lam=0.3, T0=0.5,Tf=0.01,anneal=False)
    plt.imshow(Y, cmap='gray')
    plt.title(f"Gibbs beta={Beta},Iterations={ITERA}")
    #plt.title(f"Gibbs Beta_h={Beta_h},Beta_v={Beta_v},Beta_d1={Beta_d1},Beta_d2={Beta_d2},Iterations={ITERA}")
    plt.savefig(f"Gibbs_{ITERA}")
    plt.show()
    generate_other = input('Do you want to generate another sample?(y/n)')