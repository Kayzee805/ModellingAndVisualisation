import cythonCode
import logging
import threading
import time
import numpy as np
from scipy.stats import sem
import matplotlib.pyplot as plt
size=100
sweeps=10000
pS=0.8
pI=0.1
pR=0.01

'''
Uncomment the following to run each of the task. Input args has not been provided so the user will need to type in the values for each variable manually.
'''
#cythonCode.main(size,sweeps,pS,pI,pR,runAnim=True)
#cythonCode.task3(50)
#test.task4(50,sweeps)
#cythonCode.task5(50,10000)
#cythonCode.task7(100,10000)
# #100    0.8  0.1    0.001
# #test.main(size,sweeps,pS,pI,pR,runAnim=False,genData=False,task5DO=True)

def plot5():
    array1 = np.loadtxt("data/task5/Task5_RawData1.dat")
    array2 = np.loadtxt("data/task5/Task5_RawData2.dat")
    array3 = np.loadtxt("data/task5/Task5_RawData3.dat")
    array4 = np.loadtxt("data/task5/Task5_RawData4.dat")
    array5 = np.loadtxt("data/task5/Task5_RawData5.dat")
    combined = np.array((array1,array2,array3,array4,array5))
    np.savetxt("data/Task5_RawDataCombined.dat",np.transpose(combined),fmt='%.6f')

def plot5Combined():
    allArray=np.loadtxt("data/Task5_RawDataCombined.dat")
    finalArray=[]
    errors=[]
    for i in range(len(allArray)):
        finalArray.append(np.mean(allArray[i]))
        errors.append(sem(allArray[i]))
    pImmune = np.linspace(0,1,101)
    np.savetxt("data/Task5_ProcessedData.dat",np.transpose(np.array((pImmune,finalArray,errors))),fmt='%.6f')

    plt.scatter(pImmune,finalArray,s=5,color='k')
    plt.plot(pImmune,finalArray,color='b')
    plt.errorbar(pImmune,finalArray,yerr=errors,ecolor='r')
    plt.xlabel("Immune Fraction")
    plt.ylabel("Infection Fraction")
    plt.title("Immune vs Infection fractions for 50x50 model")
    plt.savefig("figures/Task5_50x50.png")
    plt.show()


# plot5()
# plot5Combined()


def plot7():
    array1 = np.loadtxt("data/task7/Task7_RawData1.dat")
    array2 = np.loadtxt("data/task7/Task7_RawData2.dat")
    array3 = np.loadtxt("data/task7/Task7_RawData3.dat")
    array4 = np.loadtxt("data/task7/Task7_RawData4.dat")
    array5 = np.loadtxt("data/task7/Task7_RawData5.dat")
    combined = np.array((array1,array2,array3,array4,array5))
    np.savetxt("data/Task5Part2_RawDataCombined.dat",np.transpose(combined),fmt='%.6f')

def plot7Combined():
    allArray=np.loadtxt("data/Task5Part2_RawDataCombined.dat")
    finalArray=[]
    errors=[]
    for i in range(len(allArray)):
        finalArray.append(np.mean(allArray[i]))
        errors.append(sem(allArray[i]))
    pImmune = np.linspace(0,1,101)
    np.savetxt("data/Task5Part2_ProcessedData.dat",np.transpose(np.array((pImmune,finalArray,errors))),fmt='%.6f')

    plt.scatter(pImmune,finalArray,s=10,color='k')
    plt.plot(pImmune,finalArray,color='b')
    plt.errorbar(pImmune,finalArray,yerr=errors,ecolor='r')
    plt.xlabel("Immune Fraction")
    plt.ylabel("Infection Fraction")
    plt.title("Immune vs Infection fractions for 100x100 model")
    plt.savefig("figures/Task5_100x100.png")
    plt.show()
# plot7()
# plot7Combined()
def plot4():
    array1 = np.loadtxt("data/Task4_RawData1.dat")
    array2 = np.loadtxt("data/Task4_RawData2.dat")
    array3 = np.loadtxt("data/Task4_RawData3.dat")
    array4 = np.loadtxt("data/Task4_RawData4.dat")
    array5 = np.loadtxt("data/Task4_RawData5.dat")
    variance = []
    errors = []
    for i in range(len(array1)):
        temp= [array1[i],array2[i],array3[i],array4[i],array5[i]]
        variance.append(np.mean(temp))
        errors.append(sem(temp))

    # p1s = np.linspace(0.2,0.5,31)
    # print(len(variance))
    # plt.scatter(p1s,variance,s=10,color='r')
    # plt.plot(p1s,variance,color='b')
    # plt.errorbar(p1s,variance,yerr=errors,ecolor='k')
    # plt.xlabel("P1")
    # plt.ylabel("Variance")
    # plt.title("Variance vs P1 with p2=p3=0.5")
    # plt.savefig("figures/Task4_ScaledVariance_cutP3.png")
    # plt.show()
    combinedRaw = np.array((array1,array2,array3,array4,array5))
    np.savetxt("data/Task4_RawDataCombined.dat",np.transpose(combinedRaw),fmt='%.6f')

def plot4Combined():
    allArray = np.loadtxt("data/Task4_RawDataCombined.dat")
    variance =[]
    errors=[]
    for i in range(len(allArray)):
        variance.append(np.mean(allArray[i]))
        errors.append(sem(allArray[i]))
    p1s = np.linspace(0.2,0.5,31)
    print(len(variance))
    plt.scatter(p1s,variance,s=10,color='r')
    plt.plot(p1s,variance,color='b')
    plt.errorbar(p1s,variance,yerr=errors,ecolor='k')
    plt.xlabel("P1")
    plt.ylabel("Variance")
    plt.title("Variance vs P1 with p2=p3=0.5")
    plt.savefig("figures/Task4_ScaledVariance_cutP3.png")
    plt.show()
    np.savetxt("data/Task4_ProcessedData.dat",np.transpose(np.array((p1s,variance,errors))),fmt='%.6f')

#plot4()
#plot4Combined()
